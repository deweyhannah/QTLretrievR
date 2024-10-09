#' Generate mapping data and peaks for QTL analysis
#'
#' @param outdir Output directory where files mapping and peaks lists should be saved.
#' @param peaks_out String indicating the name for output peaks file. Should end in .rds.
#' @param map_out String indicating the name for output mapping file. Should end in .rds.
#' @param genoprobs Either a string with the name of the genoprobs file, or the genoprobs object.
#' @param samp_meta Sample metadata. Either a string pointing to the file, or the object itself.
#' @param expr_mats Vector of expression matrix files. One for each tissue, in the order that tissues were supplied to genoprobs.
#' @param covar_factors Additive covariate factors. These need to be columns in the factor metadata.
#' @param thrA Minimum reported LOD threshold for autosomes. Default is 5.
#' @param thrX Minimum reported LOD threshold for X chromosome. Default is 5.
#' @param gridFile File location for genome grid. Defaults to object loaded with package for 75k grid.
#' @param localRange What is defined as "local". Default is 10e6.
#' @param biomart String pointing to annotations file or annotations object.
#' @param total_cores Number of available cores to use for parallelization. Default is NULL.
#'
#' @return A list containing: \itemize{
#'  \item{maps_list}{A list of dataframes and lists to that can be used for future analyses and in other functions \itemize{
#'  \item{qtlprobs}{Genome probabilities in qtl2format}
#'  \item{covar_list}{list of covariate matrices for each tissue}
#'  \item{expr_list}{Original normalized expression data for each tissue}
#'  \item{exprZ_list}{Rank Z normalized expression data for each tissue}
#'  \item{kinship_loco}{Kisnhip Matrix calculated using the "loco" option in `qtl2::calc_kinship`}
#'  \item{gmap}{Genomic map of markers}
#'  \item{map_dat2}{Combined genomic and physical map of markers}
#'  \item{pmap}{Physical map of markers}
#'  \item{tissue_samp}{Metadata broken down for each tissue}}}
#'  \item{peaks_list}{A list of peaks list for each tissue.}}
#' @export
#'
#' @importFrom tidyr separate
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble remove_rownames column_to_rownames lst
#' @importFrom utils read.delim
#' @importFrom stats model.matrix formula
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
mapQTL <- function(outdir, peaks_out, map_out, genoprobs, samp_meta, expr_mats, covar_factors, thrA = 5, thrX = 5, gridFile = gridfile, localRange = 10e6,
                   biomart, total_cores = NULL) {
  ## Expression Matrices should be listed in the same order as tissues were for tsv2genoprobs call
  ## Load probs
  if (is.list(genoprobs)) {
    tmp_probs <- check_data(genoprobs, "genoprobs")
  }
  if (is.character(genoprobs)) {
    tmp_probs <- check_data(paste0(outdir, "/", genoprobs), "genoprobs")
  }
  probs_list <- tmp_probs
  rm(tmp_probs)

  ## Check inputs

  if (length(expr_mats) == 0) {
    stop("Please provide at least one expression matrix")
  }
  if (length(names(probs_list)) > length(expr_mats)) {
    stop("More probabilites found than expression matrices. Are you missing any expression inputs?")
  }
  if (length(names(probs_list)) < length(expr_mats)) {
    stop("More expression matrices than probabilites found. Are you using the correct probabilities?")
  }
  if (length(covar_factors) == 0) {
    stop("Please provide at least one covariate for model")
  }

  ## Modify Probs and Determine Kinship
  #if( is.null(total_cores)) total_cores <- get_cores()
  qtlprobs <- list()
  kinship_loco <- list()
  for (tissue in names(probs_list)) {
    message(tissue)
    qtlprobs[[tissue]] <- probs_3d_to_qtl2(probs_list[[tissue]])
    kinship_loco[[tissue]] <- qtl2::calc_kinship(qtlprobs[[tissue]], "loco", cores = 1)
  }

  ## Create maps
  # message(gridFile)
  if (is.character(gridFile)) {
    grid_map <- read.delim(gridFile, stringsAsFactors = F, row.names = 1)
  } else {
    grid_map <- gridFile
  }
  map_dat <- grid_map[, c("chr", "cM")]

  markers <- dimnames(probs_list[[1]])[[3]]
  map_dat <- map_dat[markers, ]
  # assert_that(are_equal(rownames(map_dat), markers))
  gmap <- split_map(map_dat)
  map_dat$marker <- rownames(map_dat)

  message("gmap complete")

  map_dat2 <- map_dat |>
    tidyr::separate(marker, into = c("chrom", "pos_bp"), convert = T, remove = F) |>
    dplyr::mutate(n = 1:dplyr::n()) |>
    tibble::as_tibble()
  pmap <- split_map(dplyr::select(
    map_dat2, marker,
    chr, pos_bp
  ) |>
    as.data.frame() |>
    tibble::remove_rownames() |>
    tibble::column_to_rownames("marker"))

  message("pmap complete")

  ## Load expression matrices
  expr_list <- list()
  for (i in 1:length(expr_mats)) {
    expr <- expr_mats[i]
    tissue <- names(probs_list)[[i]]
    expr_list[[tissue]] <- as.matrix(read.delim(expr, row.names = 1, header = T, stringsAsFactors = F))
  }

  message("expression loaded")

  ## Sample Details
  sample_details <- read.delim(samp_meta, stringsAsFactors = F, header = T)
  if (!all(covar_factors %in% colnames(sample_details))) {
    stop("Chosen factors are not in sample metadata. Please check factors and sample metadata for missing or misspelled elements")
  }
  for (fact in covar_factors) {
    sample_details[, fact] <- as.factor(sample_details[, fact])
  }

  ## Reorganize and calculate rankZ for expression matrices
  exprZ_list <- list()
  for (tissue in names(expr_list)) {
    samps_keep <- intersect(rownames(probs_list[[tissue]]), colnames(expr_list[[tissue]]))
    message(paste0("Working with ", length(samps_keep), " samples and ", nrow(expr_list[[tissue]])," genes in ", tissue, "."))
    expr_list[[tissue]] <- expr_list[[tissue]][, samps_keep, drop = FALSE]
    exprZ_list[[tissue]] <- apply(expr_list[[tissue]], 1, rankZ)
  }

  std <- standardize(expr_list, exprZ_list, sample_details, tissues = names(expr_list))
  expr_list <- std$expr
  exprZ_list <- std$exprZ
  tissue_samp <- std$tissue_samp

  message("rankZ normalized")

  ## Calculate covariate matrices
  covar_list <- list()
  for (tissue in names(tissue_samp)) {
    covar_list[[tissue]] <- model.matrix(formula(paste0("~", paste0(covar_factors, collapse = "+"))), data = tissue_samp[[tissue]])
    rownames(covar_list[[tissue]]) <- tissue_samp[[tissue]]$ID
    covar_list[[tissue]] <- covar_list[[tissue]][, colnames(covar_list[[tissue]]) != "(Intercept)", drop = FALSE]
    for (fact in covar_factors) {
      if (length(levels(tissue_samp[[tissue]][, fact])) == 2) {
        colnames(covar_list[[tissue]])[grepl(fact, colnames(covar_list[[tissue]]))] <- fact
      }
    }
  }

  message("covariates calculated")

  outfile <- paste0(outdir, "/", map_out)

  maps_list <- tibble::lst(qtlprobs, covar_list, expr_list, exprZ_list, kinship_loco, gmap, map_dat2, pmap, tissue_samp)
  # names(maps_list) <- c("qtlprobs", "covar_list", "expr_list", "exprZ_list", "kinship_loco", "gmap", "map_dat2", "pmap", "tissue_samp")
  saveRDS(maps_list, file = outfile)

  message("map saved")
  ## Run Batchmap

  message("calculating peaks")

  ## Add functionality to save out files if wanted -- probably in the background functions file.
  peaks_list <- list()

  if( is.null(total_cores)) total_cores <- get_cores()
  max_genes <- max(sapply(exprZ_list, ncol)) # Calculate the maximum number of rows across all data frames in exprZ_list
  num_tissues <-  length(names(exprZ_list))
  if( max_genes < 1000){
    cores_needed <- 8 # Limiting #of cores if there are <1000 genes in total
  }else{
    cores_needed <- total_cores
  }

  doParallel::registerDoParallel(cores = min(total_cores, cores_needed)) # no need for a lot of cores if there aren't that many genes!
  each_tissue <- floor(  min(total_cores, cores_needed) / num_tissues ) # Divide cores per tissue and pass onto the foreach loop

  message(paste0("Registering ", min(total_cores, cores_needed), " cores and passing ", each_tissue ," cores per tissue to ", num_tissues ," tissue(s). Does that look right? If not please set total_cores parameter to the number of available cores." ) )

  peak_tmp <- foreach::foreach(tissue = names(exprZ_list)) %dopar% {
    batch_wrap(
      tissue,
      exprZ_list,
      kinship_loco,
      qtlprobs,
      covar_list,
      gmap,
      thrA,
      thrX,
      each_tissue
    )
  }
  doParallel::stopImplicitCluster()

  for (i in 1:length(peak_tmp)) {
    tissue <- peak_tmp[[i]]$tissue
    peaks_list[[tissue]] <- peak_tmp[[i]]$peaks
    # message(paste0(tissue, colnames(peaks_list[[tissue]]), sep = " "))
  }

  message("adding annotations to peaks")
  peaks_list <- annotatePeaks(maps_list, peaks_list, biomart, localRange)


  outfile <- paste0(outdir, "/", peaks_out)
  saveRDS(peaks_list, file = outfile)
  message("saved peaks")

  map_peaks <- tibble::lst(maps_list, peaks_list)
  # names(map_peaks) <- c("maps_list", "peaks_list")
  return(map_peaks)
}
