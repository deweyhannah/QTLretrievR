#' Calculate founder effects for significant peaks
#'
#' @param peaks List of dataframes containing QTL peaks for each tissue
#' @param mapping List of relevant mapping information for overall project
#' @param suggLOD Suggestive LOD threshold. Default is 7.
#' @param outdir String of path to output directory where effects lists will be saved.
#' @param outfile Output file name to save mediation results for later use
#' @param n.cores Number of cores passed to qtl2. Default is 4.
#'
#' @return A list containing: \itemize{
#'  \item{effects_blup}{QTL effect BLUPs from scan along one chromosome. Output from [qtl2::scan1blup()]}
#'  \item{effects_std}{QTL effect coefficients from scan along one chromosme. Output from [qtl2::scan1coef()]}
#'  \item{peaks}{annotated peaks with LOD scores above suggestive threshold}}
#' @export
#'
#' @importFrom stats setNames
#' @importFrom dplyr mutate select rename filter arrange
#' @importFrom GenomicRanges GRanges
#' @importFrom qtl2 scan1blup scan1coef
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallel detectCores
#' @importFrom tibble lst
#'

qtl_effects <- function(mapping, peaks, suggLOD = 8, outdir, outfile, n.cores = 4) {
  ## Load in data
  if ((is.character(peaks) & is.list(mapping)) | (is.list(peaks) & is.character(mapping))) {
    stop("Peaks and mapping must both direct to an RDS file or be lists")
  }

  ## Load and organize relevant data
  if (is.list(peaks)) {
    tmp_peaks <- check_data(peaks, type = "peaks")
    tmp_map <- check_data(mapping)
  }
  if (is.character(peaks)) {
    tmp_peaks <- check_data(paste0(outdir, "/", peaks), type = "peaks")
    tmp_map <- check_data(paste0(outdir, "/", mapping))
  }
  # list2env(tmp_peaks,.GlobalEnv)
  # list2env(tmp_map,.GlobalEnvp)
  # rm(c(tmp_peaks,tmp_map))

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

  if (!is.null(tmp_peaks)) {
    peaks_list <- tmp_peaks
    rm(tmp_peaks)
  }

  message("data checked")

  peaks2 <- list()
  for (tissue in names(peaks_list)) {
    peaks2[[tissue]] <- peaks_list[[tissue]] |>
      interp_bp(df = ., genmap = gmap, physmap = pmap) |>
      dplyr::mutate(interp_bp_peak = ifelse(interp_bp_peak == 3e6, 3000001, interp_bp_peak))

    # message(paste0(names(peaks2[[tissue]]), sep = " "))
  }

  marker_list <- list()
  for (tissue in names(peaks2)) {
    marker_list[[tissue]] <- map_dat2 |>
      dplyr::select(chrom, pos_bp) |>
      dplyr::rename(start = pos_bp) |>
      dplyr::mutate(end = start) |>
      GenomicRanges::GRanges()
  }

  peaksf <- list()
  for (tissue in names(peaks2)) {
    # message(paste0(names(peaks2[[tissue]]), sep = " "))
    peaksf[[tissue]] <- peaks2[[tissue]] |>
      dplyr::filter(lod > suggLOD) |>
      dplyr::arrange(peak_chr, interp_bp_peak)
    # message(paste0(names(peaksf[[tissue]]), sep = " "))
    query <- peaksf[[tissue]] |>
      dplyr::select(peak_chr, interp_bp_peak) |>
      dplyr::rename(chrom = peak_chr, start = interp_bp_peak) |>
      dplyr::mutate(end = start) |>
      GenomicRanges::GRanges()
    subject <- marker_list[[tissue]]
    peaksf[[tissue]]$before <- map_dat2$marker[GenomicRanges::follow(query, subject)]
    peaksf[[tissue]]$after <- map_dat2$marker[GenomicRanges::precede(query, subject)]
    peaksf[[tissue]] <- peaksf[[tissue]] |>
      dplyr::filter(!is.na(before) & !is.na(after))

    # message(paste0(names(peaksf[[tissue]]), sep = " "))
    # message(paste0(head(peaksf[[tissue]]$peak_chr), sep = " "))
  }

  message("peaks extracted, calculating effects now")

  each_tissue <- floor( as.numeric(parallelly::availableCores()) / length(names(peaksf)))
  doParallel::registerDoParallel(cores = each_tissue )
  effects_out <- foreach::foreach(names(peaksf))  %dopar% {
    call_effects(
      tissue, peaksf[[tissue]], qtlprobs[[tissue]],
      gmap, exprZ_list[[tissue]], kinship_loco[[tissue]],
      covar_list[[tissue]], n.cores
    )
  }

  message("effects calculated. saving to RDS")
  # message("out names")
  # message(paste0(names(effects_out), sep = " "))
  # message("out structure")
  # message(paste0(str(effects_out), sep = " "))
  # message(head(effects_out[[1]]$effects_blup))
  # names(effects_out) <- names(peaksf)
  effects_blup <- list()
  effects_std <- list()
  for (i in 1:length(effects_out)) {
    # message(i)
    tissue <- effects_out[[i]]$tissue
    effects_blup[[tissue]] <- effects_out[[i]]$effects_blup
    effects_std[[tissue]] <- effects_out[[i]]$effects_std
  }
  peaks <- peaksf

  effects_out <- tibble::lst(effects_blup, effects_std, peaks)
  saveRDS(effects_out, paste0(outdir, "/", outfile))

  return(effects_out)
}

call_effects <- function(tissue, peaks, probs, gmap, exprZ, kinship, covars, cores) {
  message(tissue)
  n_peaks <- nrow(peaks)
  # message(n_peaks)
  # message(paste0(colnames(peaks), sep = " "))
  haps <- LETTERS[1:8]

  # message(paste0(peaks[[tissue]]$peak_chr[1:5], sep = " "))

  # message("blup")
  effects_blup <- BiocParallel::bplapply(1:n_peaks, function(i) blup_scan(i, peaks, probs, gmap, exprZ, kinship, covars, cores),
    BPPARAM = BiocParallel::MulticoreParam(workers = cores)
  )
  # message("coef")
  effects_coef <- BiocParallel::bplapply(1:n_peaks, function(i) coef_scan(i, peaks, probs, gmap, exprZ, kinship, covars)) # ,
  # BPPARAM = BiocParallel::MulticoreParam(workers = cores))

  effects_blup_tmp <- do.call("rbind", effects_blup)
  effects_coef_tmp <- do.call("rbind", effects_coef)

  effects_blup <- colMeans(effects_blup_tmp[, haps])
  effects_std <- colMeans(effects_coef_tmp[, haps])

  effect_out <- tibble::lst(tissue, effects_blup, effects_std)
  return(effect_out)
}

blup_scan <- function(i, peaks, probs, gmap, exprZ, kinship, covars, cores) {
  this_chrom <- peaks$peak_chr[i]
  # message(this_chrom)
  this_markers <- c(peaks$before[i], peaks$after[i])
  probs_2marker <- subset_probs(probs, this_chrom, this_markers)
  g <- setNames(list(gmap[[this_chrom]][this_markers]), this_chrom)
  pheno <- peaks$phenotype[i]

  out_blup <- qtl2::scan1blup(probs_2marker,
    exprZ[, pheno, drop = FALSE],
    kinship[[this_chrom]],
    covars,
    cores = cores
  )
  return(out_blup)
}

coef_scan <- function(i, peaks, probs, gmap, exprZ, kinship, covars) {
  this_chrom <- peaks$peak_chr[i]
  # message(this_chrom)
  this_markers <- c(peaks$before[i], peaks$after[i])
  probs_2marker <- subset_probs(probs, this_chrom, this_markers)
  g <- setNames(list(gmap[[this_chrom]][this_markers]), this_chrom)
  pheno <- peaks$phenotype[i]

  out_std <- qtl2::scan1coef(probs_2marker,
    exprZ[, pheno, drop = FALSE],
    kinship = kinship[[this_chrom]],
    addcovar = covars
  ) # ,
  # cores = cores)
  return(out_std)
}
