#' Plot QTL maps (peak vs gene)
#'
#' @param map_dat Mapping information for each marker used to determine genoprobs
#' @param peaks List of annotated peaks for each tissue
#' @param sigLOD Significant LOD threshold. Default is 8.
#' @param outdir String to output directory where plots should be saved
#' @param outbase File name to save plots to
#' @param psave Whether or not to save plots to png. Default is TRUE.
#' @param unit Units for start/stop/midpoint from annotations. One of "bp" or "mbp".
#' Default is "bp"
#'
#' @return A list of peak maps (ggplot objects) for each tissue
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes scale_x_discrete scale_y_discrete expansion theme element_text element_blank ggsave geom_rect geom_point scale_fill_manual
#' @importFrom ggpubr theme_pubclean
#' @importFrom tibble tibble
#' @importFrom dplyr select rename filter mutate
#'
plot_eqtlmap <- function(map_dat, peaks, sigLOD = 8, outdir, outbase, psave = T, unit = "bp"){
  if(psave == TRUE){
    if(is.null(outdir) & !is.null(outbase)){
      stop("Attempting to save plot and missing output directory")
    }
    if(!is.null(outdir) & is.null(outbase)){
      stop("Attempting to save plot and missing base file name")
    }
    if(!is.null(outdir) & !is.null(outbase)){
      message(paste0("saving plots to: ", outdir, "/", outbase, "_lod", sigLOD, "_<tissue>.png"))
    }
  }

  ## Set up chromosome midpoints and offset
  uchr <- c(as.character(1:19), "X")
  cl <- dplyr::select(map_dat, chr, pos_bp) |>
    group_by(chr) |>
    summarize(len=max(pos_bp))
  clp <- with(cl, setNames(len, chr))
  chrom_lens <- setNames(as.numeric(clp[uchr]), uchr)
  chrom_lens_offset <- cumsum(chrom_lens) - chrom_lens
  chrom_lens_midpt <- chrom_lens_offset + chrom_lens/2

  ## Set up chromosome segmentation
  chroms <-names(chrom_lens)
  chrom_segments <- tibble( start = 0,
                            end = chrom_lens,
                            chr = chroms,
                            type = as.character(rep(c(0,1),10)))
  chrom_segments$start <- chrom_segments$start+ chrom_lens_offset[chrom_segments$chr]
  chrom_segments$end <- chrom_segments$end+ chrom_lens_offset[chrom_segments$chr]

  ## Generate the plots to a list of all tissues. Option to save as png.
  for(tissue in names(peaks)){
    if(unit == "mbp"){
      peaks[[tissue]] <- peaks[[tissue]] |>
        dplyr::mutate(start = 1e6*start,
                      end = 1e6*end,
                      midpoint = 1e6*midpoint)
    }
    peaks[[tissue]] <- peaks[[tissue]] |>
      dplyr::mutate(cumsum_bp_peak = interp_bp_peak + chrom_lens_offset[peaks[[tissue]]$peak_chr],
                    cumsum_bp_gene = midpoint + chrom_lens_offset[peaks[[tissue]]$chr])
  }

  peak_map <- list()
  for(tissue in names(peaks)){
    eqtl_map <- ggplot2::ggplot()+
      ## Add vertical rectangles to distinguish between each chromosome
      ggplot2::geom_rect( data = chrom_segments, ggplot2::aes( xmin = start, xmax = end, ymin = 0, ymax = max(end), fill = type),
                 inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c("dark gray","white")) +
      ## Add peaks data and filter to significant
      ggplot2::geom_point(data = peaks[[tissue]] |>
                   dplyr::filter( lod > sigLOD),
                   ggplot2::aes( x = cumsum_bp_peak, y = cumsum_bp_gene),
                 size = 2,
                 col = "blue3",
                 inherit.aes = FALSE ) +
      ggpubr::theme_pubclean(base_size = 16) +
      ggplot2::scale_x_discrete( name = "eQTL peak",
                        limits = chrom_lens_midpt,
                        labels = names(chrom_lens),
                        expand = ggplot2::expansion( mult = 0.02))+
      ggplot2::scale_y_discrete( name = "Gene midpoint",
                                 limits = chrom_lens_midpt,
                                 labels = names(chrom_lens),
                                 expand = ggplot2::expansion( mult = 0.02)) +
      ggplot2::theme( axis.text = ggplot2::element_text(size = 10),
             panel.grid.major.x = ggplot2::element_blank(),
             panel.grid.major.y = ggplot2::element_blank())

    ## Save to list for return
    peak_map[[tissue]] <- eqtl_map

    ## Save file if wanted
    if(psave == TRUE){
      ggplot2::ggsave(paste0(outbase, "_lod", sigLOD, "_", tissue, ".png"), eqtl_map, device = png, path = outdir)
    }
  }
  return(peak_map)
}
