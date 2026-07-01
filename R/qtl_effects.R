#' Calculate founder effects for significant peaks
#'
#' @param peaks List of dataframes containing annotated QTL peaks for each
#'  tissue, or full path to `.rds` containing one.
#' @param mapping Mapping list from `mapQTL`, or full path to `.rds`
#'  containing one.
#' @param suggLOD Suggestive LOD threshold to use for filtering phenotypes.
#'  Default is 6.
#' @param outdir Directory to save effects output files. Default is `NULL`.
#' @param effects_out String indicating the name of the output file containing
#'  founder haplotype effects results. This file will be saved in `.rds` format
#'   and used for downstream analysis and visualization. Should end in `.rds`.
#'    Default is "`effects.rds`"
#' @param total_cores Number of available cores to use for parallelization.
#'  Default is `NULL`.
#' @param save Indicates object return/save behavior. One of
#'  `c("sr", "so", "ro")`; save & return, save only, return only.
#'   Default is "sr".
#' @param BPPARAM BiocParallel Parameter
#'
#' @return A list containing:
#'  \item{effects_blup}{QTL effect BLUPs from scan along one chromosome. Output
#'  from [qtl2::scan1blup()]}
#'  \item{peaks}{annotated peaks with LOD scores above suggestive threshold}
#'  These two dataframes have identical row orders, use `cbind` to merge them.
#'
#' @export
#'
#' @importFrom stats setNames
#' @importFrom dplyr mutate select rename filter arrange
#' @importFrom GenomicRanges GRanges follow precede nearest
#' @importFrom qtl2 scan1blup scan1coef
#' @importFrom BiocParallel SerialParam MulticoreParam bpnworkers
#' @importFrom tibble lst
#'

qtl_effects <- function(peaks, mapping, suggLOD = 6, outdir = NULL,
                        effects_out = "effects.rds", total_cores = NULL,
                        save = "sr", BPPARAM = BiocParallel::SerialParam()) {
  ## Load in data
  if ((is.character(peaks) & is.list(mapping)) | (is.list(peaks) &
                                                  is.character(mapping))) {
    stop("Peaks and mapping must both direct to an RDS file or be lists")
  }

  if (inherits(BPPARAM, "SerialParam") && is.null(total_cores)) {
    message("No BPPARAM provided - Detecting core availabiltiy to run parallel processes")
    n_cores <- get_cores()
    workers <- max(1, n_cores - 1)
    BPPARAM <- BiocParallel::MulticoreParam(workers = workers)
  }
  if (!is.null(total_cores)) {
    BPPARAM <- MulticoreParam(workers = total_cores)
  }

  workers <- BiocParallel::bpnworkers(BPPARAM)

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
    if ("peak_bp" %notin% colnames(peaks_list[[tissue]])) {
      peaks_list[[tissue]] <- interp_bp(df = peaks_list[[tissue]],
                                        genmap = gmap, physmap = pmap)
    }
    peaks2[[tissue]] <- peaks_list[[tissue]] |>
      dplyr::mutate(peak_bp = ifelse(peak_bp == 3e6, 3000001, peak_bp))

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
    peaksf[[tissue]] <- peaks2[[tissue]] |>
      dplyr::filter(lod >= suggLOD) |>
      dplyr::arrange(peak_chr, peak_bp)
    query <- peaksf[[tissue]] |>
      dplyr::select(peak_chr, peak_bp) |>
      dplyr::rename(chrom = peak_chr, start = peak_bp) |>
      dplyr::mutate(end = start) |>
      GenomicRanges::GRanges()
    subject <- marker_list[[tissue]]

    # peaksf[[tissue]] <- peaksf[[tissue]] |>
    #   dplyr::mutate(
    #     before = map_dat2$marker[GenomicRanges::follow(query, subject)],
    #     after  = map_dat2$marker[GenomicRanges::precede(query, subject)],
    #     marker = map_dat2$marker[GenomicRanges::nearest(query, subject)]
    #   )

    follow_idx   <- GenomicRanges::follow(query, subject)
    precede_idx  <- GenomicRanges::precede(query, subject)
    nearest_idx  <- GenomicRanges::nearest(query, subject)

    # Convert whatever comes back to a plain integer vector
    follow_idx  <- if (isS4(follow_idx))  S4Vectors::subjectHits(follow_idx)  else follow_idx
    precede_idx <- if (isS4(precede_idx)) S4Vectors::subjectHits(precede_idx) else precede_idx
    nearest_idx <- if (isS4(nearest_idx)) S4Vectors::subjectHits(nearest_idx) else nearest_idx

    peaksf[[tissue]] <- peaksf[[tissue]] |>
      dplyr::mutate(
        before = map_dat2$marker[follow_idx],
        after  = map_dat2$marker[precede_idx],
        marker = map_dat2$marker[nearest_idx]
      )

    peaksf[[tissue]] <- peaksf[[tissue]] |>
      dplyr::filter(!is.na(before) & !is.na(after))

    if(anyNA(peaksf[[tissue]]$phenotype)) {
      index <- which(is.na(peaksf[[tissue]]$phenotype))
      message(paste0(peaksf[[tissue]][index,], collapse = "\t"))
      stop("NAs found in phenotype")
    }

    peaksf[[tissue]] <- as.data.frame(peaksf[[tissue]])
  }

  message("peaks extracted, calculating effects now")

  # if( is.null(total_cores)) total_cores <- get_cores()
  # available_cores <- get_cores()
  # if( total_cores > available_cores) total_cores <- available_cores

  max_peaks <- max(vapply(peaksf, nrow, integer(1))) # how many peaks are there?
  num_tissues <-  length(names(peaksf)) # number of tissues
  cores_per_tissue <- max(1, floor(workers / num_tissues))

  effects_res <- BiocParallel::bplapply(
    names(peaksf),
    FUN = function(tissue) {
      peaks <- peaksf[[tissue]]
      call_effects(
        tissue  = tissue,
        peaks   = peaks,
        probs   = qtlprobs[[tissue]],
        gmap    = gmap,
        exprZ   = exprZ_list[[tissue]],
        kinship = kinship_loco[[tissue]],
        covars  = covar_list[[tissue]],
        cores   = cores_per_tissue
      )
    },
    BPPARAM = BPPARAM
  )

  message("effects calculated")
  effects_blup <- list()
  for (i in seq_len(length(effects_res))) {
    tissue <- effects_res[[i]]$tissue
    effects_blup[[tissue]] <- effects_res[[i]]$effects_blup
    haps <- colnames(effects_blup[[tissue]])
    # effects_blup[[tissue]] <- cbind(effects_blup[[tissue]], peaksf[[tissue]])
    # effects_blup[[tissue]] <- effects_blup[[tissue]][,c("phenotype", haps), drop = FALSE]
  }
  peaks <- peaksf

  effects_ret <- tibble::lst(effects_blup, peaks)

  if(save %in% c("sr","so")) {
    message("Saving to RDS")
    saveRDS(effects_ret, paste0(outdir, "/", effects_out))
  }
  if(save %in% c("sr","ro")) {
    return(effects_ret)
  }
}

call_effects <- function(tissue, peaks, probs, gmap, exprZ, kinship, covars,
                         cores) {
  n_peaks <- nrow(peaks)
  haps <- colnames(probs[[1]])

  ## not paralelizing scan1blup since qtl2 already does that. Instead I am
  ## passing all cores to the function to let qtl2 do the parallelization.
  if (n_peaks > 0) {
    # effects_blup <- lapply(1:n_peaks, function(i) {
    #   blup_scan(i, peaks, probs, gmap, exprZ, kinship, covars, cores = cores)
    # })
    effects_blup <- lapply(1:n_peaks, function(i) {
      result <- blup_scan(i, peaks, probs, gmap, exprZ, kinship, covars, cores = cores)
      result[peaks$marker[i], , drop = FALSE]  # take only the target marker row
    })
  } else {
    effects_blup <- data.frame(matrix(nrow = 1, ncol = length(haps)))
    colnames(effects_blup) <- haps
    effect_out <- tibble::lst(tissue, effects_blup)
    return(effect_out)
  }


  effects_blup_tmp <- do.call("rbind", effects_blup)

  effects_blup <- (effects_blup_tmp[, haps])

  effect_out <- tibble::lst(tissue, effects_blup)
  return(effect_out)
}

blup_scan <- function(i, peaks, probs, gmap, exprZ, kinship, covars, cores) {

  this_chrom <- peaks$peak_chr[i]
  this_marker <- peaks$marker[i]
  probs_2marker <- subset_probs(probs, this_chrom,this_marker )
  g <- setNames(list(gmap[[this_chrom]][this_marker]), this_chrom)
  pheno <- peaks$phenotype[i]

  if (is.null(colnames(exprZ))) {
    stop("exprZ has NULL colnames inside worker")
  }

  if (!is.matrix(exprZ)) {
    stop(sprintf("exprZ is not matrix: class = %s", paste(class(exprZ), collapse=", ")))
  }

  if (!(pheno %in% colnames(exprZ))) {
    stop(sprintf("pheno '%s' not found; ncol=%d", pheno, ncol(exprZ)))
  }


  out_blup <- qtl2::scan1blup(probs_2marker,
    exprZ[, pheno, drop = F],
    kinship[[this_chrom]],
    covars,
    cores = cores
  )
  return(out_blup)
}

