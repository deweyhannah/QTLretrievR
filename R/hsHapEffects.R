#' Map founder haplotype effects for a given hotspot.
#'
#' @param effects Effects object `qtl_effects`.
#' @param tbands List of transband locations (`bands.rna` list from
#'  `transbands` function)
#' @param chromosome Chromosome that the transband (hotspot) is present on.
#' @param tissue Tissue that the transband (hotspot) is present in
#' @param sigLOD Significant LOD threshold to use for filtering phenotypes.
#' @param hsNum If there are multiple hotspots on a chromosome,
#' indicate which one. Default is 1.
#' @param pop One of `c("do", "cc", "other")` to indicate founder population.
#'  Default is "`do`".
#' @param founders If `pop == "other"`, list of founders in haplotype order.
#' @param palette Founder color map.
#' @param topFeats Optional. Top hotspot features based on PCA analysis.
#' @param psave Logical. Save the plot as `.png`. Default `TRUE`.
#' @param pname File name to save plot (needs to end in `.png`). Default is
#'  `haplotype_effects_<tissue>_transband_<hsNum>_chromosome_<chromosome>.png`
#' @param outdir Directory to save plots. Default is `NULL`.
#'
#' @export
#'
hsHapEffects <- function(effects, tbands, chromosome, tissue, sigLOD, hsNum = 1,
                         pop = "do", founders = NULL, palette = NULL,
                         topFeats = NULL, psave = TRUE, pname = NULL,
                         outdir = NULL) {
  if (!is.null(founders) & length(founders) > 8 & is.null(palette)) {
    stop(paste0(length(founders), " founders detected,
                please provide a palette that contains at
                least that many colors"))
  }

  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (psave & !is.null(pname)) {
    pname <- paste0("haplotype_effects_", tissue, "_transband_",
                    hsNum, "_chromosome", chromosome, ".png")
    message(paste0("Plot to be saved. Saving as ", pname, " in ", outdir))
  }

  hSpot <- tbands[[tissue]] |>
    dplyr::filter(chr == chromosome)

  hsStart <- hSpot$start[hsNum]
  hsEnd <- hSpot$end[hsNum]

  if (psave & is.null(pname)) {
    pname <- paste0("haplotype_effects_", chromosome, "_", hsStart, "_",
                    hsEnd, ".png")
    message(paste0("Plot to be saved, but name not provided. Saving as ",
                   pname, " in ", outdir))
  }

  if (pop %in% c("do", "cc")) {
    founders <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
    palette <- c("#F0E442", "#555555", "#E69F00", "#0072B2",
                 "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

  }

  peak_effects <- cbind(effects$peaks[[tissue]], effects$effects_blup[[tissue]])

  peak_sub <- peak_effects[which(peak_effects$peak_chr == chromosome &
                                   peak_effects$local == 0 &
                                   peak_effects$lod > sigLOD &
                                   peak_effects$peak_bp < hsEnd &
                                   peak_effects$peak_bp > hsStart),]

  if (!is.null(topFeats)) {
    peak_sub <- peak_sub |> dplyr::filter(phenotype %in% topFeats)
  }

  rownames(peak_sub) <- NULL

  # return(peak_sub)
  num_founders <- length(founders)

  hs_mat <- peak_sub |>
    # (\(x) { rownames(x) <- NULL; x })() |>
    dplyr::select(symbol, LETTERS[seq_len(num_founders)]) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(var = "symbol") |>
    as.matrix() |>
    t()

  # return(hs_mat)

  names(palette) <- founders
  names(founders) <- LETTERS[seq_len(num_founders)]

  row_annot <- ComplexHeatmap::rowAnnotation(df = data.frame(Founders = founders),
                                             col = list(Founders = palette))

  ht <- ComplexHeatmap::Heatmap(hs_mat,
                                name = "Haplotype Effects",
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                rect_gp = grid::gpar(col = "white", lwd = 2),
                                row_names_gp = grid::gpar(fontsize = 10),
                                row_order = LETTERS[seq_len(num_founders)],
                                right_annotation = row_annot,
                                column_names_gp = grid::gpar(fontsize = 10),
                                row_title = "Founders",
                                column_title = "Target",
                                column_title_side = "top",
                                row_title_gp = grid::gpar(fontsize = 20),
                                column_title_gp = grid::gpar(fontsize = 20),
                                show_column_names = TRUE,
                                show_row_names = FALSE)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "left")

  if (psave == TRUE) {
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

}
