#' Determine significant and suggestive LOD thresholds
#'
#' @description
#' Determine significant and suggestive LOD thresholds for your QTLs, using
#' `scan1perm` from `qtl2`. Select the number of genes to permute, the number of
#'  permutations, and the number of genes to run in each batch.
#'
#'
#' @param mapping Mapping list from `mapQTL`, or full path to `.rds`
#'  containing one.
#' @param tissue Tissue to determine thresholds for.
#' @param annots Data frame of phenotype annotations for filtering to
#'  autosomal phenotypes. Columns must include "id", "symbol", "start", "end".
#' @param n.gene Number of phenotypes to run permutations on. Default is 100.
#' @param n.perm Number of permutations to run per phenotype.  Default is 500.
#' @param batch.size Number of genes in each parallelized batch.
#'
#' @return A list containing:
#' \describe{
#'   \item{perms}{A `scan1perm` object of permutations for all selected genes.}
#'   \item{thresholds}{A list with two elements:
#'     \describe{
#'       \item{significant}{Median LOD threshold at alpha = 0.05.}
#'       \item{suggestive}{Median LOD threshold at alpha = 0.4.}
#'     }
#'   }
#' }
#' @export
#'
#' @importFrom qtl2 scan1perm
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom dplyr filter select slice_sample
#' @importFrom tibble lst
#'
LOD_thld <- function(mapping, tissue, annots = NULL, n.gene = 100, n.perm = 500,
                     batch.size = 5) {
  ## Check that mapping object is valid
  if (is.list(mapping)) {
    tmp_map <- check_data(mapping)

  }
  if (is.character(mapping)) {
    tmp_map <- check_data(mapping)
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

  ## Parallelization calculations
  tot_cores <- get_cores()

  num_batches <- ceiling(n.gene / batch.size)

  cores_per_batch <- min(8, floor(tot_cores / num_batches))

  if (tot_cores <= 8 * num_batches) {
    reg_cores <- tot_cores
  } else {
    reg_cores <- 8 * num_batches
  }

  ## Prepare for permutation calculations
  if (!is.null(annots)) {
    annotations <- annots |>
      dplyr::filter(id %in% colnames(exprZ_list[[tissue]]),
                    chr %in% c(1:19)) |>
      dplyr::select(id) |>
      dplyr::slice_sample(n = n.gene)
    genes <- annotations$id
  } else {
    genes <- colnames(exprZ_list[[tissue]])
    genes <- sample(genes, size = n.gene)
  }
  batch_genes <- split(genes, ceiling(seq_along(genes) / batch.size))

  message(paste0("Using ", n.gene, " phenotypes: ", paste0(genes,
                                                           collapse = "  ")))

  message(paste0("Running ", n.perm, " permutations on ", n.gene, " genes in ",
                 num_batches, " batches."))

  ## Run permutation calculations
  doParallel::registerDoParallel(cores = reg_cores)
  perms_out2 <- foreach::foreach(gene = batch_genes) %dopar% {
    qtl2::scan1perm(genoprobs = qtlprobs[[tissue]],
                    pheno     = exprZ_list[[tissue]][, gene],
                    kinship   = kinship_loco[[tissue]],
                    addcovar  = covar_list[[tissue]],
                    n_perm    = n.perm,
                    cores     = cores_per_batch)
  }
  doParallel::stopImplicitCluster()

  ## Prepare data for return
  perms_df <- do.call(cbind, perms_out2)

  sigLOD <- stats::median(summary(perms_df, alpha = 0.05))
  suggLOD <- stats::median(summary(perms_df, alpha = 0.4))

  thlds <- tibble::lst(sigLOD, suggLOD)
  out <- tibble::lst(perms_df, thlds)

  return(out)
}
