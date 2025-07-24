#' Generate mapping data and peaks for QTL analysis
#' @description
#' This function performs QTL mapping across multiple tissues, generating genome-wide LOD scores and identifying significant peaks. It supports parallel processing and flexible input formats for expression data and covariates.
#'
#'
#' @param genoprobs A `qtl2`-formatted genotype probabilities object in list form, one probabilities object per tissue, or a character path to a saved `.rds` file containing one.
#' @param samp_meta Sample metadata. Either a string pointing to the file, or the object itself.
#' @param expr_mats List of expression matrices (objects), or character paths to the file. One matrix per tissue. The order *must match* the tissue order in `genoprobs`
#' @param covar_factors Character vector of column names in `samp_meta` to be used as additive covariates.
#' @param thrA Minimum reported LOD threshold for autosomes. Default is 5.
#' @param thrX Minimum reported LOD threshold for X chromosome. Default is 5.
#' @param gridFile File location for genome grid. Defaults to object loaded with package for 75k grid.
#' @param localRange What is defined as "local". Default is 10e6.
#' @param outdir Output directory where files mapping and peaks lists should be saved. Default is NULL.
#' @param peaks_out String indicating the name for output peaks file. Should end in ".rds". Default is "peaks.rds"
#' @param map_out String indicating the name for output mapping file. Should end in ".rds". Default is "map.rds"
#' @param annots String pointing to annotations file or annotations object.
#' @param total_cores Number of available cores to use for parallelization. Default is NULL.
#' @param save Character. Determines output behavior: `"sr"` to save and return, `"so"` to save only, `"ro"` to return only. Default is `"sr"`.
#' @param rz Logical. Set to `TRUE` if expression data is already rank Z-transformed. Default is `FALSE`.
#' @param phys Logical. if `TRUE`, use the physical map for peak calling; otherwise use the genomic map. Default is `TRUE`.
#'
#'
#' @return A list containing: \itemize{
#'  \item{maps_list}{A list of dataframes and lists to that can be used for future analyses and in other functions \itemize{
#'  \item{qtlprobs}{Genome probabilities in qtl2format}
#'  \item{covar_list}{list of covariate matrices for each tissue}
#'  \item{expr_list}{Original normalized expression data for each tissue}
#'  \item{exprZ_list}{Rank Z-transformed expression data for each tissue}
#'  \item{kinship_loco}{Kinship matrix calculated using the "loco" option in `qtl2::calc_kinship`}
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
mapQTL <- function(genoprobs, samp_meta, expr_mats, covar_factors, thrA = 5, thrX = 5, gridFile = gridfile,
                   localRange = 10e6, outdir = NULL, peaks_out = "peaks.rds", map_out = "map.rds",
                   annots = NULL, total_cores = NULL, save = "sr", rz = F, phys = T) {

  ## Check save conflicts
  if (save %in% c("sr", "so")) {
    if (is.null(outdir)) {
      stop("Requested Save. No output directory provided, and no default.")
    }
  }

  ## Expression Matrices should be listed in the same order as tissues in genoprobs list
  ## Load probs
  if (is.list(genoprobs)) {
    tmp_probs <- check_data(genoprobs, type = "genoprobs")
  }
  if (is.character(genoprobs)) {
    tmp_probs <- check_data(genoprobs, type = "genoprobs")
  }
  probs_list <- tmp_probs
  rm(tmp_probs)

  if (!is.null(annots)) {
    if (is.character(annots)) {
      annots <- read.delim(annots, stringsAsFactors = F, header = T)
    }
  }

  ## Check inputs
  ## Assumes one expression matrix per tissue, matching the order of tissues in genoprobs

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

  ## Convert genotype probabilities and calculate LOCO kinship matrices for each tissue
  qtlprobs <- probs_list
  kinship_loco <- list()
  for (tissue in names(probs_list)) {
    message(tissue)
    kinship_loco[[tissue]] <- qtl2::calc_kinship(qtlprobs[[tissue]], "loco", cores = 1)
  }

  ## Create maps
  if (is.character(gridFile)) {
    grid_map <- read.delim(gridFile, stringsAsFactors = F, row.names = 1)
  } else {
    grid_map <- gridFile
  }
  map_dat <- grid_map[, c("chr", "cM", "pos")]

  markers <- c()
  for (chrom in names(qtlprobs[[1]])) {
    markers <- c(markers, dimnames(qtlprobs[[1]][[chrom]])[[3]])
  }
  map_dat <- map_dat[markers, ]
  gmap <- split_map(map_dat)
  map_dat$marker <- rownames(map_dat)

  message("gmap complete")

  map_dat2 <- map_dat |>
    dplyr::mutate(n = 1:dplyr::n(), chrom = chr) |>
    dplyr::rename(pos_bp = pos) |>
    tibble::as_tibble()

  pmap <- split_map(dplyr::select(
    map_dat2, marker,
    chr, pos_bp
  ) |>
    as.data.frame() |>
    tibble::remove_rownames() |>
    tibble::column_to_rownames("marker"))

  message("pmap complete")

  ## Load expression matrices for each tissue. If input is a file path, rad it; otherwise assume its already a matrix object
  expr_list <- list()
  for (i in 1:length(expr_mats)) {
    expr <- expr_mats[i]
    tissue <- names(probs_list)[[i]]
    if (is.character(expr)) {
      expr_list[[tissue]] <- as.matrix(read.delim(expr, row.names = 1, header = T, stringsAsFactors = F))
    }
    else {
      expr_list[[tissue]] <- expr[[tissue]]
    }
  }

  if (rz) {
    message("rankZ transformed counts provided")
    for (tissue in names(expr_list)) {
      expr_list[[tissue]] <- t(expr_list[[tissue]])
    }
  }

  message("expression loaded")

  if (!is.null(annots)) {
    for (tissue in names(expr_list)) {
      if (length(intersect(rownames(expr_list[[tissue]]), annots$id)) == 0) {
        stop(paste0("Phenotypes in ", tissue, " not present in annotations file. Please check phenotype names in counts."))
      }
      if (length(intersect(rownames(expr_list[[tissue]]), annots$id)) < nrow(expr_list[[tissue]])) {
        message(paste0("Not all phenotypes in ", tissue, " are present in annotations file. Not all phenotypes will be annotated."))
      }
    }
  }

  ## Sample Details
  if (is.character(samp_meta)) {
    sample_details <- read.delim(samp_meta, stringsAsFactors = F, header = T)
  }
  if (is.data.frame(samp_meta)) {
    sample_details <- samp_meta
  }
  if (!all(covar_factors %in% colnames(sample_details))) {
    stop("Chosen factors are not in sample metadata. Please check factors and sample metadata for missing or misspelled elements")
  }
  for (fact in covar_factors) {
    sample_details[, fact] <- as.factor(sample_details[[fact]])
  }

  ## Normalize expression data using rankZ transformation (unless already transformed). Align samples with genotype data.
  exprZ_list <- list()
  if (rz) {
    for (tissue in names(expr_list)) {
      exprZ_list[[tissue]] <- t(expr_list[[tissue]])
    }
  } else {
    for (tissue in names(expr_list)) {
      samps_keep <- intersect(rownames(probs_list[[tissue]][[1]]), colnames(expr_list[[tissue]]))
      message(paste0("Working with ", length(samps_keep), " samples and ", nrow(expr_list[[tissue]])," genes in ", tissue, "."))
      expr_list[[tissue]] <- expr_list[[tissue]][, samps_keep, drop = FALSE]
      exprZ_list[[tissue]] <- apply(expr_list[[tissue]], 1, rankZ)
    }
  }

  std <- standardize(expr_list, sample_details, tissues = names(expr_list))
  expr_list <- std$expr
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

  maps_list <- tibble::lst(qtlprobs, covar_list, expr_list, exprZ_list, kinship_loco, gmap, map_dat2, pmap, tissue_samp)

  if(save %in% c("sr","so")) {
    outfile <- paste0(outdir, "/", map_out)
    saveRDS(maps_list, file = outfile)
    message("map saved")
  }

  ## Run Batchmap
  message("calculating peaks")

  peaks_list <- list()

  available_cores <- get_cores()
  if( is.null(total_cores)) total_cores <- available_cores
  if( total_cores > available_cores) total_cores <- available_cores
  max_genes <- max(sapply(exprZ_list, ncol)) # Calculate the maximum number of rows across all data frames in exprZ_list
  num_tissues <-  length(names(exprZ_list))
  if( max_genes < 1000){
    cores_needed <- 8 # Limiting # of cores if there are <1000 genes in total
  }else{
    cores_needed <- total_cores
  }

  if (phys) {
    map <- pmap
  } else {
    map <- gmap
  }

  doParallel::registerDoParallel(cores = min(total_cores, cores_needed)) # no need for a lot of cores if there aren't that many genes!

  ## Divide available cores evenly across tissues for parallel peak calling
  each_tissue <- floor(  min(total_cores, cores_needed) / num_tissues )

  message(paste0("Registering ", min(total_cores, cores_needed), " cores and passing ", each_tissue ," cores per tissue to ", num_tissues ," tissue(s). Does that look right? If not please set total_cores parameter to the number of available cores." ) )

  peak_tmp <- foreach::foreach(tissue = names(exprZ_list)) %dopar% {
    batch_wrap(
      tissue,
      exprZ_list,
      kinship_loco,
      qtlprobs,
      covar_list,
      map,
      thrA,
      thrX,
      each_tissue,
      phys
    )
  }
  doParallel::stopImplicitCluster()

  for (i in 1:length(peak_tmp)) {
    tissue <- peak_tmp[[i]]$tissue
    peaks_list[[tissue]] <- peak_tmp[[i]]$peaks
  }

  if (!is.null(annots)) {
    message("adding annotations to peaks")
    if (!phys) {
      for(tissue in names(peaks_list)) {
        peaks_list[[tissue]] <- interp_bp(peaks_list[[tissue]], gmap, pmap)
      }
    }
    peaks_list <- annotatePeaks(maps_list, peaks_list, annots, localRange)
  }

  if(save %in% c("sr","so")) {
    outfile <- paste0(outdir, "/", peaks_out)
    saveRDS(peaks_list, file = outfile)
    message("saved peaks")
  }

  ## Return both mapping and peak results as a named list
  if(save %in% c("sr","ro")) {
    map_peaks <- tibble::lst(maps_list, peaks_list)
    return(map_peaks)
  }
}
