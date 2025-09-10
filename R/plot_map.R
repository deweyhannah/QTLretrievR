#' Plot QTL maps (peak vs gene)
#'
#' @param map_dat `map_dat2` from `mapQTL` mapping list.
#' @param peaks List of dataframes containing QTL peaks for each tissue
#'  (annotated).
#' @param sigLOD Significant LOD threshold to use for filtering phenotypes.
#'  Default is 7.5
#' @param psave Logical. Save the plot as `.png`. Default `TRUE`.
#' @param pname File name to save plot (needs to end in `.png`).
#'  Default is "eqtl_map_LOD<sigLOD>_<tissue>.png".
#' @param outdir Directory to save plots. Default is `NULL`.
#' @param unit One of `c("bp", "mbp")`; annotation position units.
#'  Default is "bp".
#' @param map_col Plot color. Default is "#0073C2FF"
#'
#' @return A list of peak maps (ggplot objects) for each tissue
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes scale_x_discrete scale_y_discrete expansion
#' theme element_text element_blank ggsave geom_rect geom_point
#' scale_fill_manual
#' @importFrom ggpubr theme_pubclean
#' @importFrom tibble tibble
#' @importFrom dplyr select rename filter mutate
#'
plot_eqtlmap <- function(map_dat, peaks, sigLOD = 7.5, outdir = NULL,
                         pname = NULL, psave = TRUE, unit = "bp",
                         map_col = "#0073C2FF") {
  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (psave & !is.null(pname)) {
    message(paste0("Plot to be saved. Saving as ", pname, " in ", outdir))
  }
  if (psave & is.null(pname)) {
    temp_name <- paste0("eqtl_map_LOD",sigLOD,"_<tissue>.png")
    message(paste0("Plot to be saved, but name not provided. Saving as ",
                   temp_name, " in ", outdir))

  }

  ## Set up chromosome midpoints and offset
  uchr <- c(as.character(1:19), "X")
  cl <- dplyr::select(map_dat, chr, pos_bp) |>
    group_by(chr) |>
    summarize(len = max(pos_bp))
  clp <- with(cl, setNames(len, chr))
  chrom_lens <- setNames(as.numeric(clp[uchr]), uchr)
  chrom_lens_offset <- cumsum(chrom_lens) - chrom_lens
  chrom_lens_midpt <- chrom_lens_offset + chrom_lens / 2

  ## Set up chromosome segmentation
  chroms <- names(chrom_lens)
  chrom_segments <- tibble(
    start = 0,
    end = chrom_lens,
    chr = chroms,
    type = as.character(rep(c(0, 1), 10))
  )
  chrom_segments$start <- chrom_segments$start +
    chrom_lens_offset[chrom_segments$chr]
  chrom_segments$end <- chrom_segments$end +
    chrom_lens_offset[chrom_segments$chr]

  ## Generate the plots to a list of all tissues. Option to save as png.
  for (tissue in names(peaks)) {
    if (unit == "mbp") {
      peaks[[tissue]] <- peaks[[tissue]] |>
        dplyr::mutate(
          start = 1e6 * start,
          end = 1e6 * end,
          midpoint = 1e6 * midpoint
        )
    }
    peaks[[tissue]] <- peaks[[tissue]] |>
      dplyr::mutate(
        cumsum_bp_peak = peak_bp + chrom_lens_offset[peaks[[tissue]]$peak_chr],
        cumsum_bp_gene = midpoint + chrom_lens_offset[peaks[[tissue]]$chr]
      )
  }

  peak_map <- list()
  for (tissue in names(peaks)) {
    eqtl_map <- ggplot2::ggplot() +
      ## Add vertical rectangles to distinguish between each chromosome
      ggplot2::geom_rect(
        data = chrom_segments, ggplot2::aes(xmin = start, xmax = end, ymin = 0,
                                            ymax = max(end), fill = type),
        inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = c("dark gray", "white")) +
      ## Add peaks data and filter to significant
      ggplot2::geom_point(
        data = peaks[[tissue]] |>
          dplyr::filter(lod > sigLOD),
        ggplot2::aes(x = cumsum_bp_peak, y = cumsum_bp_gene),
        size = 2,
        col = map_col,
        inherit.aes = FALSE
      ) +
      ggpubr::theme_pubclean(base_size = 16) +
      ggplot2::scale_x_continuous(
        name = "eQTL peak",
        breaks = chrom_lens_midpt,
        labels = names(chrom_lens),
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::scale_y_continuous(
        name = "Gene midpoint",
        breaks = chrom_lens_midpt,
        labels = names(chrom_lens),
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 16),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
      )

    ## Save to list for return
    peak_map[[tissue]] <- eqtl_map

    ## Save file if wanted
    if (psave) {
      pname <- paste0("eqtl_map_LOD",sigLOD,"_",tissue,".png")
      ggplot2::ggsave(pname, eqtl_map, device = "png", path = outdir,
                      width = 3072, height = 3072, units = "px")
    }
  }
  return(peak_map)
}
