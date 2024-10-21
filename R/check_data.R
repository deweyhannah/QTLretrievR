#' Checking the inputted data for QTLretrievR
#'
#' @param x can be a path to an Rds object or a list of R objects.
#' @param type The name of the object if a single one is passed. Use "genoprobs" for genotype probabilities, "mediation" for mediation results and "peaks" for a table of eQTL peaks.
#' @return List of properly names R objects to be used in downstream analysis of QTLretrievR.
#'
#'
#' @importFrom tibble lst
#'

check_data <- function(x, type = "") {
  ## Load data into environment

  if (is.character(x)) {
    y <- readRDS(x)
    # message(names(y))
    # message(length(intersect(names(y), c("probs_list","map_dat2","peaks_list","res_list","effects_blup"))))
    if (length(intersect(names(y), c("probs_list", "map_dat2", "peaks_list", "res_list", "effects_blup"))) == 1) {
      exprZ_list <- y$exprZ_list
      covar_list <- y$covar_list
      expr_list <- y$expr_list
      gmap <- y$gmap
      kinship_loco <- y$kinship_loco
      map_dat2 <- y$map_dat2
      pmap <- y$pmap
      qtlprobs <- y$qtlprobs
      tissue_samp <- y$tissue_samp
    }
    if (length(intersect(names(y), c("probs_list", "map_dat2", "peaks_list", "res_list", "effects_blup"))) == 0) {
      if (type == "genoprobs") {
        probs_list <- y
      }
      if (type == "mediation") {
        res_list <- y
      }
      if (type == "peaks") {
        peaks_list <- y
      }
    }
  }
  if (is.list(x)) {
    # list2env(x,.GlobalEnv)
    if (length(intersect(names(x), c("probs_list", "map_dat2", "peaks_list", "res_list", "effects_blup"))) == 1) {
      exprZ_list <- x$exprZ_list
      covar_list <- x$covar_list
      expr_list <- x$expr_list
      gmap <- x$gmap
      kinship_loco <- x$kinship_loco
      map_dat2 <- x$map_dat2
      pmap <- x$pmap
      qtlprobs <- x$qtlprobs
      tissue_samp <- x$tissue_samp
    }
    if (length(intersect(names(x), c("probs_list", "map_dat2", "peaks_list", "res_list", "effects_blup"))) == 0) {
      if (type == "genoprobs") {
        probs_list <- x
      }
      if (type == "mediation") {
        res_list <- x
      }
      if (type == "peaks") {
        peaks_list <- x
      }
    }
  }

  ## Run checks
  if (type == "genoprobs") {
    stopifnot(is.list(probs_list))
    for (tissue in names(probs_list)) {
      stopifnot(names(probs_list[[tissue]]) == c(1:19,"X"))
    }
    return(probs_list)
  }
  if (type = "") {
    ## Class check
    stopifnot(is.list(qtlprobs))
    stopifnot(is.list(covar_list))
    stopifnot(is.list(expr_list))
    stopifnot(is.list(exprZ_list))
    stopifnot(is.list(kinship_loco))
    stopifnot(is.list(gmap))
    stopifnot(is.list(pmap))
    stopifnot(is.list(tissue_samp))
    stopifnot(is.data.frame(map_dat2))

    ## Names (Tissue) check
    a <- data.frame(table(c(
      names(qtlprobs), names(covar_list), names(expr_list),
      names(exprZ_list), names(kinship_loco), names(tissue_samp)
    )))
    stopifnot(length(unique(a$Freq)) == 1)

    ## Names (chromosome) check
    b <- data.frame(table(c(names(kinship_loco[[1]]), names(gmap), names(pmap), names(qtlprobs[[1]]))))
    stopifnot(length(unique(b$Freq)) == 1)

    ## Sample Check
    for (tissue in names(tissue_samp)) {
      stopifnot(identical(sort(tissue_samp[[tissue]]$ID), sort(colnames(expr_list[[tissue]]))))
      stopifnot(identical(sort(tissue_samp[[tissue]]$ID), sort(rownames(exprZ_list[[tissue]]))))
      stopifnot(identical(sort(tissue_samp[[tissue]]$ID), sort(rownames(covar_list[[tissue]]))))
      stopifnot(identical(sort(rownames(expr_list[[tissue]])), sort(colnames(exprZ_list[[tissue]]))))
    }

    ## map_dat2 column check
    stopifnot(all(c("chr", "cM", "marker", "chrom", "pos_bp", "n") %in% colnames(map_dat2)))

    return(tibble::lst(qtlprobs, covar_list, expr_list, exprZ_list, kinship_loco, gmap, pmap, tissue_samp, map_dat2))
  }
  if (type = "peaks") {
    stopifnot(is.list(peaks_list))
    for (tissue in names(peaks_list)) {
      stopifnot(is.data.frame(peaks_list[[tissue]]))
      stopifnot(all(c("phenotype", "peak_chr", "peak_cM", "lod", "ci_lo", "ci_hi") %in% colnames(peaks_list[[tissue]])))
    }
    return(peaks_list)
  }
  if (type = "mediation") {
    stopifnot(is.list(res_list))
    for (tissue in names(res_list)) {
      stopifnot(all(c("target_id", "qtl_lod", "qtl_chr", "mediator", "mediator_id", "mediator_chr", "mediator_midpoint", "LOD") %in% colnames(res_list[[tissue]])))
    }
    return(res_list)
  }
  if (type = "effects") {
    stopifnot(is.list(effects_blup))
    stopifnot(is.list(effects_std))

    stopifnot(identical(sort(names(effects_blup)), sort(names(effects_blup))) &
      identical(sort(names(effects_blup)), sort(names(peaks))))

    for (tissue in names(effects_blup)) {
      stopifnot(is.matrix(effects_blup[[tissue]]))
      stopifnot(is.matrix(effects_std[[tissue]]))

      stopifnot(colnames(effects_blup[[tissue]]) == LETTERS[1:8])
      stopifnot(colnames(effects_std[[tissue]]) == LETTERS[1:8])
    }
  }
  return(tibble::lst(effects_blup, effects_std, peaks))
}
