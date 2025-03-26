#' Plot peaks associated with specific genes.
#'
#' @param mapping Mapping list from `mapQTL`
#' @param tissue Tissue to derive plots from
#' @param pheno Phenotype to derive plots from
#' @param saveDir Directory to save plots to
#' @param annots Annotations (example `annot_105`), required if plotting effects on chromosome
#' @param effects Boolean - Plotting effects on chromosome, default FALSE
#'
#' @return None
#' @export
#'
#' @importFrom qtl2 scan1 scan1blup plot_coefCC
#'
peak_plot <- function(mapping, tissue, pheno, saveDir, annots = NULL, effects = FALSE) {
  if (is.list(mapping)) {
    tmp_map <- check_data(mapping)

  }
  if (is.character(mapping)) {
    tmp_map <- check_data(paste0(outdir, "/", mapping))
  }

  if (!is.null(tmp_map)) {
    exprZ_list <- tmp_map$exprZ_list
    covar_list <- tmp_map$covar_list
    expr_list <- tmp_map$expr_list
    gmap <- tmp_map$gmap
    kinship_loco <- tmp_map$kinship_loco
    map_dat2 <- tmp_map$map_dat2
    pmap <- tmp_map$pmap
    qtlprobs <- tmp_map$qtlprobs
    tissue_samp <- tmp_map$tissue_samp
  }
  rm(tmp_map)
  rm(mapping)

  if (is.null(annots) & effects) {
    stop("Annotations needed to produce effects plot")
  }

  avail_cores <- get_cores()

  scan1_out <- qtl2::scan1(genoprobs = qtlprobs[[tissue]],
                           pheno     = exprZ_list[[tissue]][, pheno, drop = FALSE],
                           kinship   = kinship_loco[[tissue]],
                           addcovar  = covar_list[[tissue]],
                           cores     = min(4, avail_cores))

  if (effects) {
    chr <- annots$chr[annots$id == pheno]
    c2eff <- qtl2::scan1blup(genoprobs = qtlprobs[[tissue]][,as.character(chr)],
                             pheno     = exprZ_list[[tissue]][, pheno, drop = FALSE])

    png(paste0(saveDir, "/", tissue, "_", pheno, "_effects.png"), width = 1024, height = 1024)
    qtl2::plot_coefCC(x            = c2eff,
                      map          = pmap,
                      scan1_output = scan1_out,
                      bgcolor      = "gray95")
    dev.off()

    legend_plot <- paste0(saveDir,"/CC_effects_legend.png")
    if (!file.exists(legend_plot)) {
      png(legend_plot, bg = "transparent")
      plot.new()
      legend("center", c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB"), lty = 1, col = qtl2::CCcolors, bty = "n", ncol = 4, lwd = 1.3)
      dev.off()
    }

  } else {
    png(paste0(saveDir, "/", tissue, "_", pheno, "_LOD.png"), width = 1024, height = 730)
    plot(scan1_out, pmap)
    dev.off()
  }
}
