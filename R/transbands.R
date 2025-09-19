#' Identify distal hotspots above a suggestive and significant LOD score.
#'
#' @param map_dat `map_dat2` from `mapQTL` mapping list.
#' @param peaks List of dataframes containing QTL peaks for each
#'  tissue (annotated).
#' @param sigLOD Significant LOD threshold to use for filtering phenotypes.
#'  Default is 7.5
#' @param suggLOD Suggestive LOD threshold to use for filtering phenotypes.
#' Default is 6.
#' @param psave Logical. Save the plot as `.png`. Default `TRUE`.
#' @param pname File name to save plot (needs to end in `.png`).
#'  Default is `hotspots_<tissue>.png`
#' @param outdir Directory to save plots. Default is `NULL`.
#' @param color Plot color. Default is "#0073C2FF".
#'
#' @return A list containing: \itemize{
#'  \item{bands.rna}{A list of tibbles containing hotspot information for
#'   each tissue
#'  Information includes the following for each hotspot:
#'  - Chromosome of the hotspot
#'  - Start and Stop (bp) of the hotspot
#'  - The number of significant peaks within the hotspot
#'  - The number of suggestive peaks within the hotspot }
#'  \item{trans_band_plots}{A list of ggplot2 objects (one for each tissue),
#'  showing the number of significant distal hotspot peaks.}}
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_x_continuous expansion
#' xlab ylab theme element_text ggsave
#' @importFrom ggpubr theme_pubclean
#' @importFrom tibble lst as_tibble
#' @importFrom dplyr select rename group_by summarize filter mutate starts_with
#' @importFrom stats quantile
#' @importFrom GenomicRanges GRanges countOverlaps seqnames start end
#' slidingWindows nearest
#'
transbands <- function(map_dat, peaks, sigLOD = 7.5, suggLOD = 6, psave = TRUE,
                       pname = NULL, outdir = NULL, color = "#0073C2FF") {
  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (psave & is.null(pname)) {
    ptemp <- "hotspots_<tissue>.png"
    message(paste0("Plot to be saved. Saving as ", ptemp, " in ", outdir))
  }

  pname_og <- pname

  ## Set up chromosome midpoints and offset
  uchr <- c(as.character(1:19), "X")
  cl <- dplyr::select(map_dat, chr, pos_bp) |>
    group_by(chr) |>
    summarize(len = max(pos_bp))
  clp <- with(cl, setNames(len, chr))
  chrom_lens <- setNames(as.numeric(clp[uchr]), uchr)
  chrom_lens_offset <- cumsum(chrom_lens) - chrom_lens
  chrom_lens_midpt <- chrom_lens_offset + chrom_lens / 2

  ## Prep map for peaks
  ## Make chromosomes factors
  map_dat$chromF <- factor(map_dat$chrom, levels = c(as.character(1:19), "X"))
  map_dat <- map_dat |> dplyr::rename(pos_cM = cM)

  ## Organise markers into different GenomicRanges objects
  chrom_markers <- select(map_dat, chromF, n) |>
    dplyr::rename(chrom = chromF) |>
    group_by(chrom) |>
    summarize(start = min(n), end = max(n)) |>
    GenomicRanges::GRanges()

  markers_bynum <- select(map_dat, chrom, n) |>
    dplyr::rename(start = n) |>
    mutate(end = start) |>
    GenomicRanges::GRanges()

  markers <- select(map_dat, chrom, pos_bp) |>
    dplyr::rename(start = pos_bp) |>
    mutate(end = start) |>
    GenomicRanges::GRanges()

  ## Sliding windows of chromosome based markers - 50 <> windows 10 <> step
  windows <- unlist(GenomicRanges::slidingWindows(chrom_markers, width = 50,
                                                  step = 10))

  ## Start per tissue analysis with annotated peaks (suggestive only)

  message("identifying hotspots")
  bands.rna <- list()
  for (tissue in names(peaks)) {
    ## Significant distant peaks
    distant_rna <- peaks[[tissue]] |>
      dplyr::filter(lod > sigLOD, local == 0) |>
      dplyr::select(peak_chr, peak_bp) |>
      dplyr::rename(chrom = peak_chr, end = peak_bp) |>
      dplyr::mutate(start = end) |>
      GenomicRanges::GRanges()

    ## Suggestive distant peaks
    distant_rna_sugg <- peaks[[tissue]] |>
      dplyr::filter(lod > suggLOD, local == 0) |>
      dplyr::select(peak_chr, peak_bp) |>
      dplyr::rename(chrom = peak_chr, end = peak_bp) |>
      dplyr::mutate(start = end) |>
      GenomicRanges::GRanges()

    ## Identify hotspots
    hotspot <- GenomicRanges::nearest(distant_rna, markers)
    hotspot_sugg <- GenomicRanges::nearest(distant_rna_sugg, markers)

    ## Add hotspots to windows
    windows$distant_rna <- GenomicRanges::countOverlaps(windows,
                                                        markers_bynum[hotspot])
    windows$distant_rna_sugg <- GenomicRanges::countOverlaps(
      windows, markers_bynum[hotspot_sugg])

    ## Fill out data for all windows
    window_counts <- tibble(
      chrom = as.character(GenomicRanges::seqnames(windows)),
      start = GenomicRanges::start(windows),
      end = GenomicRanges::end(windows),
      distant_rna = windows$distant_rna,
      distant_rna_sugg = windows$distant_rna_sugg
    )

    mm <- match(window_counts$start, map_dat$n)
    m2 <- match(window_counts$end, map_dat$n)
    window_counts$pos_cM_start <- map_dat$pos_cM[mm]
    window_counts$pos_bp_start <- map_dat$pos_bp[mm]
    window_counts$pos_cM_end <- map_dat$pos_cM[m2]
    window_counts$pos_bp_end <- map_dat$pos_bp[m2]
    # message(paste0(colnames(window_counts), sep = " "))
    window_counts <- window_counts |>
      dplyr::mutate(midpoint = (pos_cM_end + pos_cM_start) / 2, 4)

    x <- select(window_counts, chrom, starts_with("pos_bp"),
                starts_with("distant")) |>
      filter(
        distant_rna >= quantile(distant_rna, 0.995)
      )

    ## Collapse overlapping windows into one big window
    bands <- x |>
      dplyr::rename(start = pos_bp_start, end = pos_bp_end) |>
      GenomicRanges::GRanges() |>
      GenomicRanges::reduce()

    bands$distant_rna <- GenomicRanges::countOverlaps(bands, distant_rna)
    bands$distant_rna_sugg <- GenomicRanges::countOverlaps(bands,
                                                           distant_rna_sugg)

    ## Convert to tibble and save to list
    bands.rna[[tissue]] <- bands |>
      as_tibble() |>
      mutate(chr = seqnames)
  }

  ## Make plots and save if wanted
  message("generating plots")
  trans_band_plot <- list()
  for (tissue in names(peaks)) {
    eqtl_counts <- bands.rna[[tissue]] |>
      tibble::as_tibble() |>
      dplyr::select(chrom = seqnames, start, end, distant_rna) |>
      dplyr::mutate(chrom = factor(chrom, levels = c(seq(1:19), "X"))) |>
      # adding all the marker locations to match axes
      dplyr::mutate(hotspot_midpoint = (start + end) / 2) |>
      # adding all the markers with 0 hotspot values to match axes
      rbind((map_dat |>
        dplyr::select(chrom, start = pos_bp, end = pos_bp) |>
        dplyr::mutate(
          distant_rna = 0,
          hotspot_midpoint = start
        ) |>
        dplyr::mutate(chrom = factor(chrom, levels = c(seq(1:19), "X")))))

    message(paste0(eqtl_counts$hotspot_midpoint[which(
      eqtl_counts$distant_rna != 0)], sep = " "))

    eqtl_counts$midpoint_offset <- eqtl_counts$hotspot_midpoint +
      chrom_lens_offset[eqtl_counts$chrom]

    message(paste0(unique(eqtl_counts$chrom[which(
      eqtl_counts$distant_rna != 0)]), sep = " "))
    message(paste0(eqtl_counts$midpoint_offset[which(
      eqtl_counts$distant_rna != 0)], sep = " "))

    trans_band_plot[[tissue]] <- eqtl_counts |>
      ggplot2::ggplot() +
      ggplot2::aes(x = midpoint_offset, y = distant_rna) +
      ggplot2::geom_bar(stat = "identity", width = 50, col = color,
                        fill = color) +
      ggpubr::theme_pubclean(base_size = 16) +
      ggplot2::scale_x_continuous(
        name = "Chr",
        breaks = chrom_lens_midpt,
        labels = names(chrom_lens),
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("# of distant eQTL") +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 10))

    if (psave == TRUE) {
      if (is.null(pname_og)) {
        pname <- paste0("hotspots_", tissue, ".png")
      }
      ggplot2::ggsave(pname, plot = trans_band_plot[[tissue]], device = "png",
                      path = outdir)
    }
  }
  table_plot <- tibble::lst(bands.rna, trans_band_plot)
  return(table_plot)
}
