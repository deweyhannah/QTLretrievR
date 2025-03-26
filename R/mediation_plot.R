#' Plot the change in mediation scores for the top 5
#'
#' @param effects Effects object from `qtl_effects()`
#' @param mediation Mediation object from `run_mediate()`
#' @param trans.bands Transbands object from `transbands()`
#' @param annots Molecular annotations. Default is ensembl v105 annotations.
#' @param tissue The tissue your transband of interest is in
#' @param peak.chr The chromosome your transband of interest is on
#' @param peak.num The transband number. If there is more than one transband on your chromosome of interest, specify which one (left to right) if there is only one then set this to 1
#' @param outdir Output directory where the plot will be saved. Only required if `psave = TRUE`
#' @param sigLOD Significant LOD threshold, determined by `LOD_thld()`
#' @param local How are you defining local in the transband. Default is 5e6
#' @param psave Should the plot be saved, or just returned? Default is FALSE
#'
#' @return A list containing: \itemize{
#'  \item{results}{The mediation results from your hotspot with ranks}
#'  \item{meds_ranked_sum}{A summary of the mediators and the number of times they are the top mediator for any given target}
#'  \item{p}{Mediation Plot showing the LOD drop of the top 5 mediators for each target eQTL}}
#' @export
#'
#' @importFrom dplyr inner_join filter left_join select group_by mutate arrange summarize ungroup n
#' @import ggplot2
#' @importFrom scales rescale
#' @importFrom forcats as_factor
#' @importFrom ggpubr theme_pubclean
#' @importFrom tibble lst
#'
#'
mediation_plot <- function(effects, mediation, trans.bands, annots = annot_105, tissue, peak.chr,
                     peak.num, outdir = NULL, sigLOD = 7.5, local = 5e6, psave = FALSE) {

  if (psave & is.null(outdir)) {
    stop("Please provide an output directory if you want to save the plot")
  }
  ## Merge peaks and effects into a single dataframe
  peak_effects <- cbind(effects$peaks[[tissue]], effects$effects_blup[[tissue]])

  ## Pull the peaks in the transband of interest
  band_peaks <- peak_effects |>
    dplyr::filter(peak_chr == peak.chr,
                  lod > sigLOD,
                  local == 0,
                  interp_bp_peak >= as.numeric(trans.bands$bands.rna[[tissue]][which(trans.bands$bands.rna[[tissue]]$chr == peak.chr), "start"][peak.num]),
                  interp_bp_peak <= as.numeric(trans.bands$bands.rna[[tissue]][which(trans.bands$bands.rna[[tissue]]$chr == peak.chr), "end"][peak.num]))

  ## Pull the mediation results for the peaks in the transband of interest
  peak_meds <- mediation[[tissue]] |>
    dplyr::inner_join(band_peaks |>
                        dplyr::select(target_id = phenotype,
                                      qtl_chr = peak_chr,
                                      target_pos = interp_bp_peak,
                                      target_chr = chr))

  ## Filter mediators
  meds_filt <- peak_meds |>
    dplyr::left_join(annots |> dplyr::select(target = symbol,
                                      target_id = id)) |>
    dplyr::select(target,
                  qtl_chr,
                  qtl_lod,
                  target_pos,
                  target_chr,
                  mediator_id,
                  mediator,
                  mediator_chr,
                  mediator_midpoint,
                  mediation_lod = LOD) |>
    dplyr::group_by(target, qtl_chr, qtl_lod, target_chr) |>
    dplyr::mutate(scaled_LOD = scale(mediation_lod)) |>
    dplyr::filter(abs(target_pos - mediator_midpoint) <= local,
                  qtl_chr == mediator_chr)

  ## Rank mediators
  meds_ranked <- meds_filt |>
    dplyr::mutate(mediation_lod = ifelse(target == mediator, NA, mediation_lod),
                  lod_drop = qtl_lod - mediation_lod) |>
    dplyr::group_by(target) |>
    dplyr::arrange(lod_drop) |>
    dplyr::mutate(rank = rep(seq(1:dplyr::n())))

  meds_ranked_sum <- meds_ranked |>
    dplyr::filter(rank %in% c( 1)) %>%
    dplyr::group_by(mediator) %>%
    dplyr::summarize(n = length(target), min_drop = min(lod_drop, na.rm = T), max_drop = max(lod_drop, na.rm = T), med_drop = median(lod_drop, na.rm = T)) %>%
    dplyr::arrange(desc(n))

  ## Reorganize ranked mediators for plots
  results <- meds_ranked |>
    dplyr::ungroup() |>
    dplyr::select(mediator, target, lod_drop, mediator_midpoint, target_pos, rank) |>
    dplyr::filter(rank %in% c(1:5)) |>
    dplyr::mutate(lod_drop = ifelse(lod_drop < 0, 0, lod_drop), lod_drop = ifelse(lod_drop > sigLOD, sigLOD, lod_drop)) |>
    dplyr::arrange(mediator_midpoint, mediator) |>
    dplyr::mutate(xlabel = forcats::as_factor(sort(mediator))) |>
    dplyr::arrange(target_pos) |>
    dplyr::mutate(ylabel = forcats::as_factor(target))

  ## Create plot to variable
  p <- results |>
    ggplot(aes(x = xlabel, y = ylabel)) +
    geom_point(aes(color = lod_drop, size = exp(lod_drop / 3))) +
    scale_color_gradientn(colors = c("white", "firebrick3", "navy"),
                          values = scales::rescale(c(0, 4, 8)),
                          name = "LOD\ndifference",
                          limits = c(0, 8)) +
    scale_size(breaks = 0:6,
               labels = as.character(0:6),
               range = c(0,5)) +
    guides(size = "none") +
    ggpubr::theme_pubclean(base_size = 18) +
    theme(axis.text.y = element_text(size = 8, hjust = 1),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    ylab("Target eQTL")+
    xlab("Mediator")

  ## Save the plot if psave == TRUE
  if (psave) {
    plot_title <- paste0("mediation_res_chr_", peak.chr, "_hotspot_", peak.num, ".png")
    ggsave(plot_title, device = "png", path = outdir, width = 10, height = 20, units = "in")
  }

  return(tibble::lst(results |> dplyr::select(-xlabel, -ylabel), meds_ranked_sum, p))
}
