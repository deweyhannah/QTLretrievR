#' Identify top mediators within a hotspot (molecular QTL)
#'
#' @param peaks Dataframe of peaks to pull targets from (output of `mapQTL`, then select `peaks_list$tissue`).
#' @param meds Dataframe of mediation to pull potential mediators from (output from `modiFinder`, then select `$tissue`).
#' @param tbands Dataframe of all transbands (output of `transbands`, then select `bands.rna$tissue`).
#' @param chromosome The chromosome that the hotspot falls on
#' @param hsNum If there is more than one hotspot on the chromosome, which one do you want? Default 1.
#' @param top_n The number of top mediators per target to show. Default 5.
#' @param psave Should the plot be saved? Default TRUE.
#' @param pname String indicating the name of plot to save as a .png. Default "mediation_plot_chr<chromosome>_<start>_<stop>_top_<top_n>.png".
#' @param outdir String indicating the location to save the plot to. Default NULL, required if `psave == TRUE`.
#' @param plot One of c("padj", "pval", "per_drop", "ranks") depending on what statistic should be plotted in the heatmap. Default "padj".
#'
#' @return A list containing a dataframe representing the ranking of each mediator within the hotspot for each target of the hotspot and the heatmap object.
#'
#' @importFrom dplyr filter group_by mutate ungroup arrange select between anti_join join_by case_when if_else
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom viridis viridis
#' @importFrom tidyr pivot_wider
#' @importFrom tibble lst column_to_rownames
#' @importFrom grid gpar
#' @importFrom grDevices png dev.off
#'
#' @export
#'
medPlot_hotSpot <- function(peaks, meds, tbands, chromosome, hsNum = 1, top_n = 5, psave = T, pname = NULL, outdir = NULL, plot = "padj") {
  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (psave & !is.null(pname)) {
    message(paste0("Plot to be saved. Saving as ", pname, " in ", outdir))
  }

  hSpot <- tbands |>
    dplyr::filter(chr == chromosome)

  hsStart <- hSpot$start[hsNum]
  hsEnd <- hSpot$end[hsNum]

  if (psave & is.null(pname)) {
    pname <- paste0("mediation_plot_chr", chromosome, "_", hsStart, "_", hsEnd, "_top_", top_n, ".png")
    message(paste0("Plot to be saved, but name not provided. Saving as ", pname, " in ", outdir))
  }


  feats <- peaks$phenotype[which(peaks$peak_chr == chromosome & peaks$local == 0 & peaks$peak_bp < hsEnd & peaks$peak_bp > hsStart)]

  peak_med <- hs_sig(df    = meds,
                     start = hsStart,
                     stop  = hsEnd,
                     feat  = feats,
                     chr   = chromosome)

  l2p_wide <- peak_med |>
    dplyr::group_by(target_id) |>
    dplyr::mutate(ranks = rank(-LOD_drop, ties.method = "average")) |>
    dplyr::ungroup() |>
    dplyr::mutate(padj = dplyr::case_when(ranks <= top_n ~ padj,
                                          ranks > top_n ~ NA)) |>
    dplyr::mutate(padj = dplyr::if_else(padj == 0.0000, 1e-15, padj)) |>
    dplyr::arrange(mediator_midpoint) |>
    dplyr::select(target, mediator, as.name(plot)) |>
    tidyr::pivot_wider(names_from = mediator, values_from = as.name(plot)) |>
    tibble::column_to_rownames(var = "target") |>
    as.matrix()

  ht <- ComplexHeatmap::Heatmap(l2p_wide,
                          col = rev(viridis::viridis(100)),
                          na_col = "gray95",
                          cluster_rows = F,
                          cluster_columns = F,
                          rect_gp = grid::gpar(col = "white", lwd = 2),
                          row_names_gp = grid::gpar(fontsize = 10),
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_title = "Target",
                          column_title = "Mediator",
                          column_title_side = "bottom",
                          row_title_gp = grid::gpar(fontsize = 20),
                          column_title_gp = grid::gpar(fontsize = 20),
                          heatmap_legend_param = list(title = plot),
                          show_row_names = T,
                          show_column_names = T)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "left")

  rank_return <- peak_med |>
    dplyr::group_by(target_id) |>
    dplyr::mutate(per_drop = LOD_drop / qtl_lod,
                  ranks = rank(-LOD_drop, ties.method = "average")) |>
    dplyr::ungroup()

  if (psave) {
    n_rows <- nrow(l2p_wide)
    n_cols <- ncol(l2p_wide)

    # Define base size and scaling
    base_width <- 450
    base_height <- 400
    col_scale <- 20
    row_scale <- 15

    width <- base_width + col_scale * n_cols
    height <- base_height + row_scale * n_rows

    png(paste0(outdir, "/", pname), width = width, height = height)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "left")
    dev.off()
  }

  return(tibble::lst("meds_ranked" = rank_return, "heatmap" = ht))

}


hs_filter <- function(df, start, stop, feat) {
  ## filter mediation results to +/- x Mb from peak location
  filtered_df <- df |>
    dplyr::filter(qtl_chr == mediator_chr) |>
    dplyr::filter(target_id %in% feat) |>
    dplyr::filter(dplyr::between(mediator_midpoint, left = start, right = stop))

  return(filtered_df)
}


hs_sig <- function(df, start, stop, feat, chr) {
  filtered_df <- df |> dplyr::filter(target_id %in% feat, qtl_chr == chr)
  filtered_df$LOD_drop <- filtered_df$qtl_lod - filtered_df$LOD
  peak_sp_df <- hs_filter(filtered_df, start, stop, feat)
  # peak_sp_df <- peak_sp_df |> dplyr::filter(mediator_chr == chr,
  #                                           qtl_chr == chr)

  peak_sp_df$pval <- NA

  for(feature in unique(peak_sp_df$target_id)) {
    feature_df <- filtered_df |> dplyr::filter(target_id == feature)
    feature_peak_df <- peak_sp_df |> dplyr::filter(target_id == feature)

    nullDist <- dplyr::anti_join(feature_df, feature_peak_df, by = dplyr::join_by(target_id, qtl_lod, qtl_chr, mediator, mediator_id, mediator_chr, mediator_midpoint, LOD, target, LOD_drop))
    nullDist_drop <- nullDist$LOD_drop
    for(mediator in feature_peak_df$mediator) {
      med_drop <- feature_peak_df$LOD_drop[which(feature_peak_df$mediator == mediator & feature_peak_df$qtl_chr == chr)]
      all_drops <- c(nullDist_drop, med_drop)
      ranks <- rank(all_drops, ties.method = "average")
      percentile <- ranks[length(ranks)] / length(ranks)
      z <- qnorm(percentile)
      p <- 1 - stats::pnorm(abs(z))
      peak_sp_df$pval[which(peak_sp_df$target_id == feature & peak_sp_df$mediator == mediator & peak_sp_df$qtl_chr == chr)] <- p
    }
  }
  peak_sp_df <- peak_sp_df |>
    dplyr::group_by(target_id) |>
    dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"))
  return(peak_sp_df)
}


