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
#' @param vert Logical. Rotate plot to vertical orientation. Default `FALSE`
#' @param ... Additional arguments to pass to ComplexHeatmap
#'
#'
#' @return List of plots. `ht_plot`: Plot only, formatted with legend on the left.
#' `ht`: Original ComplexHeatmap object.
#' @export
#'
hsHapEffects <- function(effects, tbands, chromosome, tissue, sigLOD, hsNum = 1,
                         pop = "do", founders = NULL, palette = NULL,
                         topFeats = NULL, psave = TRUE, pname = NULL,
                         outdir = NULL, vert = FALSE, ...) {
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
    names(founders) <- LETTERS[1:8]
    palette <- c("#F0E442", "#555555", "#E69F00", "#0072B2",
                 "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

  }

  palette_to_use <- c("#4575B4", "#4979B6", "#4E7DB8", "#5282BB", "#5786BD",
                      "#5C8BBF", "#608FC2", "#6594C4", "#6998C6", "#6E9DC9",
                      "#73A1CB", "#77A6CE", "#7CAAD0", "#80AFD2", "#85B3D5",
                      "#8AB8D7", "#8EBCD9", "#93C0DB", "#98C3DD", "#9CC6DF",
                      "#A1CAE1", "#A6CDE2", "#ABD0E4", "#B0D3E6", "#B4D6E8",
                      "#B9D9E9", "#BEDCEB", "#C3E0ED", "#C8E3EF", "#CCE6F0",
                      "#D1E9F2", "#D6ECF4", "#DBEFF6", "#E0F3F7", "#E1F3F4",
                      "#E3F4F1", "#E5F5ED", "#E7F5EA", "#E9F6E6", "#EBF7E3",
                      "#EDF8DF", "#EFF8DC", "#F0F9D8", "#F2FAD5", "#F4FBD2",
                      "#F6FBCE", "#F8FCCB", "#FAFDC7", "#FCFDC4", "#FEFEC0",
                      "#FEFEBD", "#FEFCBA", "#FEFAB7", "#FEF8B5", "#FEF6B2",
                      "#FEF4AF", "#FEF2AC", "#FEF0A9", "#FEEFA6", "#FEEDA3",
                      "#FEEBA1", "#FEE99E", "#FEE79B", "#FEE598", "#FEE395",
                      "#FEE192", "#FEDF8F", "#FDDA8C", "#FDD589", "#FDD085",
                      "#FDCB82", "#FDC67F", "#FDC17B", "#FDBC78", "#FDB775",
                      "#FCB271", "#FCAD6E", "#FCA86B", "#FCA367", "#FC9E64",
                      "#FC9961", "#FC945D", "#FC8F5A", "#FA8A57", "#F88454",
                      "#F67E51", "#F4794E", "#F1734B", "#EF6D48", "#ED6845",
                      "#EB6242", "#E85D3F", "#E6573C", "#E45139", "#E24C36",
                      "#DF4633", "#DD4030", "#DB3B2D", "#D9352A", "#D73027")

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

  founders <- factor(founders, levels = founders)


  row_annot <- ComplexHeatmap::rowAnnotation(
    df = data.frame(Founders = founders),
    col = list(Founders = setNames(palette, founders)),
    annotation_legend_param = list(
      at = founders,
      labels = founders
    )
  )



  if (!vert) {
    ht <- ComplexHeatmap::Heatmap(hs_mat,
                                  col = palette_to_use,
                                  name = "Haplotype Effects",
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
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
  }
  if (vert) {
    col_annot <- ComplexHeatmap::columnAnnotation(df = data.frame(Founders = founders),
                                               col = list(Founders = palette))
    hs_mat <- t(hs_mat)
    ht <- ComplexHeatmap::Heatmap(hs_mat,
                                  col = palette_to_use,
                                  name = "Haplotype Effects",
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  rect_gp = grid::gpar(col = "white", lwd = 2),
                                  row_names_gp = grid::gpar(fontsize = 10),
                                  column_order = LETTERS[seq_len(num_founders)],
                                  top_annotation = col_annot,
                                  column_names_gp = grid::gpar(fontsize = 10),
                                  column_title = "Founders",
                                  row_title = "Target",
                                  column_title_side = "top",
                                  row_title_gp = grid::gpar(fontsize = 20),
                                  column_title_gp = grid::gpar(fontsize = 20),
                                  show_column_names = FALSE,
                                  show_row_names = TRUE)
  }

  pdf(NULL)

  ht_plot <- ComplexHeatmap::draw(ht, heatmap_legend_side = "left",
                                  annotation_legend_side = "left",
                                  merge_legends = TRUE,
                                  legend_gap = grid::unit(5, "mm"))
  dev.off()

  if (psave == TRUE) {
    n_rows <- nrow(hs_mat)
    n_cols <- ncol(hs_mat)

    # Define base size and scaling
    base_width <- 450
    base_height <- 400
    col_scale <- 20
    row_scale <- 15
    if (vert) {
      base_width <- 400
      base_height <- 450
    }

    width <- base_width + col_scale * n_cols
    height <- base_height + row_scale * n_rows

    png(paste0(outdir, "/", pname), width = width, height = height)
    ht_plot
    dev.off()
  }

  return(tibble::lst(ht_plot, ht))
}
