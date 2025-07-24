#' Generate mapping data and peaks for G x E QTL analysis
#'
#' @description
#' This function performs G x E QTL mapping across multiple tissues, generating genome-wide LOD scores and identifying significant peaks. It supports parallel processing and flexible input formats for expression data and covariates.
#'
#'
#' @param genoprobs Either a string with the name of the genoprobs file, or the genoprobs object.
#' @param samp_meta Sample metadata. Either a string pointing to the file, or the object itself.
#' @param expr_mats Vector of expression matrix files. One for each tissue, in the order that tissues were supplied to genoprobs.
#' @param covar_factors Additive covariate factors. These need to be columns in the factor metadata.
#' @param thrA Minimum reported LOD threshold for autosomes. Default is 5.
#' @param thrX Minimum reported LOD threshold for X chromosome. Default is 5.
#' @param gridFile File location for genome grid. Defaults to object loaded with package for 75k grid.
#' @param localRange What is defined as "local". Default is 10e6.
#' @param outdir Output directory where files mapping and peaks lists should be saved. Default is NULL.
#' @param peaks_out String indicating the name for output peaks file. Should end in ".rds". Default is "gxe_peaks.rds"
#' @param map_out String indicating the name for output mapping file. Should end in ".rds". Default is "gxe_map.rds"
#' @param annots String pointing to annotations file or annotations object.
#' @param total_cores Number of available cores to use for parallelization. Default is NULL.
#' @param save Should files be saved, returned, or both. Default is "sr" (save and return). To save only use "so", to return only use "ro".
#' @param ctrl String indicating your control or background gene name (ex: "ctrl" or "CTRL")
#' @param env String indicating your exposed/treated samples (ex: "trt" or "treated" or "<your_treatment_here>")
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
#'
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
gxeQTL <- function(genoprobs, samp_meta, expr_mats, covar_factors, thrA = 5, thrX = 5, gridFile = gridfile,
                   localRange = 10e6, outdir = NULL, peaks_out = "gxe_peaks.rds", map_out = "gxe_map.rds",
                   annots = NULL, total_cores = NULL, save = "sr", delta = FALSE, ctrl, env) {

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
    tmp_probs <- check_data(paste0(outdir, "/", genoprobs), type = "genoprobs")
  }
  probs_list <- tmp_probs
  rm(tmp_probs)

  ## Check inputs
  ## Assumes one expression matrix per tissue, matching the order of tissues in genoprobs

  if (length(expr_mats) == 0) {
    stop("Please provide at least one expression matrix")
  }
  if (env %notin% names(genoprobs)) {
    stop("Please provide genoprobs for your environmental condition")
  }
  if (env %notin% names(expr_mats) & ctrl %notin% names(expr_mats)) {
    stop("Missing counts for Gene and Environment")
  }
  if (env %notin% names(expr_mats)) {
    stop("Missing counts for Environment")
  }
  if (ctrl %notin% names(expr_mats)) {
    stop("Missing counts for Gene")
  }
  if (length(covar_factors) == 0) {
    stop("Please provide at least one covariate for model")
  }

  ## Convert genotype probabilities and calculate LOCO kinship matrices for each tissue
  qtlprobs <- probs_list
  kinship_loco <- list()
  kinship_loco[[env]] <- qtl2::calc_kinship(qtlprobs[[env]], "loco", cores = 1)

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
  for (i in names(expr_mats)) {
    expr <- expr_mats[[i]]
    tissue <- i
    if (is.character(expr)) {
      expr_list[[tissue]] <- as.matrix(read.delim(expr, row.names = 1, header = T, stringsAsFactors = F))
    }
    else {
      message(paste0("passing dataframe with: ", ncol(expr), " samples, and ", nrow(expr), " genes"))
      expr_list[[tissue]] <- expr
    }
  }

  message("expression loaded")

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
    sample_details[, fact] <- as.factor(sample_details[, fact])
  }

  ## Normalize expression data using rankZ transformation (unless already transformed). Align samples with genotype data.
  exprZ_list <- list()
  for (tissue in names(expr_list)) {
    samps_keep <- intersect(rownames(probs_list[[tissue]][[1]]), colnames(expr_list[[tissue]]))
    message(paste0("Working with ", length(samps_keep), " samples and ", nrow(expr_list[[tissue]])," genes in ", tissue, "."))
    expr_list[[tissue]] <- expr_list[[tissue]][, samps_keep, drop = FALSE]
    exprZ_list[[tissue]] <- apply(expr_list[[tissue]], 1, rankZ)
  }

  std <- standardize(expr_list, exprZ_list, sample_details, tissues = names(expr_list))
  expr_list <- std$expr
  exprZ_list <- std$exprZ
  tissue_samp <- std$tissue_samp

  message("rankZ normalized")

  exprZ_gxe_genes <- intersect(colnames(exprZ_list[[env]]), colnames(exprZ_list[[ctrl]]))
  exprZ_gxe <- list()
  exprZ_gxe[[env]] <- exprZ_list[[env]][,exprZ_gxe_genes, drop = F]
  exprZ_gxe[[ctrl]] <- exprZ_list[[ctrl]][,exprZ_gxe_genes, drop = F]

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

  maps_list <- tibble::lst(qtlprobs, covar_list, expr_list, exprZ_gxe, kinship_loco, gmap, map_dat2, pmap, tissue_samp)

  if(save %in% c("sr","so")) {
    outfile <- paste0(outdir, "/", map_out)
    saveRDS(maps_list, file = outfile)
    message("map saved")
  }

  ## Run Batchmap
  message("calculating peaks")

  peaks_list <- list()

  available_cores <- get_cores()
  if( total_cores > available_cores) total_cores <- available_cores
  max_genes <- ncol(exprZ_gxe[[env]]) # Calculate the maximum number of rows across all data frames in exprZ_list
  if( max_genes < 1000){
    cores_needed <- 8 # Limiting # of cores if there are <1000 genes in total
  }else{
    cores_needed <- total_cores
  }
  #
  doParallel::registerDoParallel(cores = min(total_cores, cores_needed)) # no need for a lot of cores if there aren't that many genes!

  message(paste0(names(exprZ_gxe), collapse = "\t"))

  if (delta) {
    exprZ_delta <- list()
    exprZ_delta[[env]] <- exprZ_gxe[[env]] - exprZ_gxe[[ctrl]]
    peaks_list[[env]] <- batch_wrap(env,
                                    exprZ_delta,
                                    kinship_loco,
                                    qtlprobs,
                                    covar_list,
                                    gmap,
                                    thrA,
                                    thrX,
                                    min(total_cores, cores_needed))
  }
  else {
    peaks_list[[env]] <- batch_gxe(
      exprZ_list    = exprZ_gxe,
      kinship_loco  = kinship_loco,
      qtlprobs      = qtlprobs,
      covars        = covar_list,
      gmap          = gmap,
      thrA          = thrA,
      thrX          = thrX,
      cores         = min(total_cores, cores_needed),
      ctrl          = ctrl,
      env           = env,
      covar_factors = covar_factors
    )
  }

  doParallel::stopImplicitCluster()

  if (!is.null(annots)) {
    message("adding annotations to peaks")
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

  return(peaks_list)
}

peak_gxe <- function(i,ss, exprZ, kinship_loco, genoprobs, covar, gmap, thrA = 5, thrX = 5, n.cores = 4, ctrl, env, covar_factors) {
  start <- ss[i] + 1
  end <- ss[i + 1]

  exprZ_env <- exprZ[[env]][,start:end, drop = F]
  exprZ_ctrl <- exprZ[[ctrl]][,colnames(exprZ_env), drop = F]

  out <- lapply(1:ncol(exprZ_env), function(x){
    if(x %% 100 == 0){message(paste0("  --Trait ",x," out of ",ncol(exprZ_env)))}
    ctrl_exp <- exprZ_ctrl[,x]
    # gather covariates
    covar <- cbind(covar[[env]], ctrl_exp)
    gxe_scan1_out <- qtl2::scan1(genoprobs = genoprobs,
                                 pheno     = exprZ_env[,x],
                                 kinship   = kinship_loco,
                                 addcovar  = covar,
                                 cores     = n.cores)
    colnames(gxe_scan1_out) <- dimnames(exprZ_env)[[2]][x]
    return(gxe_scan1_out)
  })
  out <- Reduce(cbind,out)

  peaks <- qtl2::find_peaks(out, gmap,
                            drop = 1.5,
                            threshold = thrA, thresholdX = thrX
  )
  peaks <- peaks %>%
    dplyr::select(-lodindex) %>%
    dplyr::rename(phenotype = lodcolumn, peak_chr = chr, peak_cM = pos)

  return(peaks)
}

batch_gxe <- function(exprZ_list, kinship_loco, qtlprobs,
                       covars, gmap, thrA, thrX, cores, ctrl, env, covar_factors) {

  num.batches <- max(c(round(ncol(exprZ_list[[env]])/1000), 2))
  nn <- ncol(exprZ_list[[env]])
  ss <- round(seq(0, nn, length.out=num.batches + 1))

  message(paste0("Mapping ", nn, " genes in ", num.batches + 1, " batches"))

  # Calculate the maximum number of concurrent batches
  # Each batch should at least have 4 cores.
  max_concurrent_batches <- max(1, floor(cores / 4) )
  # if the #of batches > max_concurrent_batches adjust the cores
  if( num.batches > max_concurrent_batches){
    # Adjust the number of cores to use based on concurrency limit
    cores_to_use <- min(cores, max_concurrent_batches * 4)
    # get the cores to use per batch, minimum 4
    cores_per_batch <- max(4, floor(cores_to_use/max_concurrent_batches))
  } else{
    # can use all the cores
    cores_to_use <- cores
    # get the cores to use per batch
    cores_per_batch <- floor(cores_to_use/num.batches)
  }

  doParallel::registerDoParallel(cores = cores_to_use)

  # Initialize an empty list to store the results
  all_results <- list()
  for (i in seq(1, num.batches, by = max_concurrent_batches)) {
    # Determine the current batch range
    current_batches <- i:min(i + max_concurrent_batches - 1, num.batches)

    # Run the current batch of tasks in parallel
    current_results <- foreach::foreach( j = current_batches)  %dopar% {
      peak_gxe(i             = j,
               ss            = ss,
               exprZ         = exprZ_list,
               kinship_loco  = kinship_loco[[env]],
               genoprobs     = qtlprobs[[env]],
               covar         = covars,
               gmap          = gmap,
               thrA          = thrA,
               thrX          = thrX,
               n.cores       = cores_per_batch,
               ctrl          = ctrl,
               env           = env,
               covar_factors = covar_factors)
    }

    # Append the current results to the overall results
    all_results <- c(all_results, current_results)
  }

  doParallel::stopImplicitCluster()

  peaks <- do.call("rbind", all_results)

  return(peaks)
}
