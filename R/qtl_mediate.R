#' Prepare and run mediation for a set of QTL peaks
#'
#' @param peaks List of dataframes containing QTL peaks for each tissue
#' @param mapping List of relevant mapping information for overall project
#' @param suggLOD Suggestive LOD to use as filter for mediation. Default is 7.
#' @param outdir String of path to output directory where mediation lists will be saved.
#' @param annots String pointing to annotations file or annotations object.
#' @param med_out Output file name to save mediation results for later use
#' @param total_cores Number of available cores to use for parallelization. Default is NULL.
#' @param save Should files be saved, returned, or both. Default is "sr" (save and return). To save only use "so", to return only use "ro".
#'
#' @return A list containing mediation results for each tissue
#'
#' @export
#'
#' @importFrom dplyr rename filter mutate select
#' @importFrom purrr compact
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
run_mediate <- function(peaks, mapping, suggLOD = 7, outdir, annots, med_out, total_cores = NULL, save = "sr") {
  # mediate_env <- new.env()
  message("load annotations")
  if (is.character(annots)) {
    annots <- read.delim(annots)
  }
  if ("Gene.start..bp." %in% colnames(annots)) {
    message("renaming annots columns")
    annots <- annots |>
      dplyr::rename(id = Gene.stable.ID, symbol = MGI.symbol, start = Gene.start..bp., end = Gene.end..bp., chr = Chromosome.scaffold.name)
  }
  if ("gene.id" %in% colnames(annots)) {
    colnames(annots)[which(colnames(annots) == "gene.id")] <- "id"
  }
  if ("end" %in% colnames(annots)) {
    colnames(annots)[which(colnames(annots) == "end")] <- "stop"
  }
  if ("chr" %in% colnames(annots)) {
    colnames(annots)[which(colnames(annots) == "chr")] <- "chromosome"
  }

  message("checking peaks and mapping")
  if ((is.character(peaks) & is.list(mapping)) | (is.list(peaks) & is.character(mapping))) {
    stop("Peaks and mapping must both direct to an RDS file or be lists")
  }

  ## Load and organize relevant data
  if (is.list(peaks)) {
    peaks_list <- check_data(peaks, type = "peaks")
    tmp_map <- check_data(mapping)
    # if(!is.null(tmp_map)){
    #   list2env(mapping, .GlobalEnv)
    # }
  }
  if (is.character(peaks)) {
    peaks_list <- check_data(paste0(outdir, "/", peaks), type = "peaks")
    tmp_map <- check_data(paste0(outdir, "/", mapping))
    # message(paste0(names(tmp_map), sep = " "))
  }

  #stopifnot(str(peaks_list[[1]]) == "data.frame")
  if (!is.data.frame(peaks_list[[1]])) {
    stop("The elements of 'peaks_list' must be data frames.")
  }

  message("data checked")

  ## clean up
  # list2env(tmp_peaks,.GlobalEnv)
  # list2env(tmp_map,.GlobalEnv)
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

  annots_list <- list()
  for (tissue in names(peaks_list)) {
    annots_list[[tissue]] <- annots |>
      dplyr::filter(id %in% rownames(expr_list[[tissue]])) |>
      dplyr::mutate(midpoint = (start + stop) / 2)
  }

  # message(paste0(ls(), sep = " "))

  ## Identify tissues and targets
  qtl_peaks <- list()
  qtl_target <- list()
  for (tissue in names(peaks_list)) {
    qtl_peaks[[tissue]] <- peaks_list[[tissue]] |>
      dplyr::filter(lod > suggLOD) |>
      # interp_bp(.) |>
      dplyr::mutate(phenotype = gsub("_.*", "", phenotype)) |>
      dplyr::mutate(target_id = phenotype) |>
      dplyr::filter(target_id %in% annots_list[[tissue]]$id)
    qtl_target[[tissue]] <- exprZ_list[[tissue]][, annots_list[[tissue]]$id]
  }

  message("filtered peaks")

  ## Run annotations
  targ_annot <- list()
  for (tissue in names(peaks_list)) {
    # message(paste0(names(annots_list[[tissue]]), sep = " "))
    targ_annot[[tissue]] <- annots_list[[tissue]] |>
      dplyr::mutate(target_id = id, chrom = chromosome) |>
      dplyr::mutate(chrom = ifelse(chromosome == "MT", "M", chromosome)) |>
      dplyr::filter(!is.na(chrom)) |>
      dplyr::mutate(chr = chrom, pos = abs(stop + start) / 2)
  }
  targ_covar <- covar_list
  probs <- qtlprobs
  kinship <- kinship_loco

  qtl_mediatior <- list()
  med_annot <- list()
  for (tissue in names(exprZ_list)) {
    med_annot[[tissue]] <- targ_annot[[tissue]] |>
      dplyr::rename(mediator_id = target_id)
    qtl_mediatior[[tissue]] <- exprZ_list[[tissue]][, annots_list[[tissue]]$id]
    qtl_mediatior[[tissue]] <- qtl_mediatior[[tissue]][, med_annot[[tissue]]$id, drop = FALSE]
    qtl_target[[tissue]] <- qtl_target[[tissue]][, targ_annot[[tissue]]$target_id, drop = FALSE]
    qtl_peaks[[tissue]] <- qtl_peaks[[tissue]] |>
      dplyr::filter(target_id %in% targ_annot[[tissue]]$target_id)
  }
  med_covar <- covar_list

  message("running mediation")

  if( is.null(total_cores)) total_cores <- get_cores()
  available_cores <- get_cores()
  if( total_cores > available_cores) total_cores <- available_cores
  max_peaks <- max(sapply(qtl_peaks, nrow)) # get the maximum number of peaks
  num_tissues <-  length(names(qtl_peaks)) # number of tissues
  if( max_peaks < 1000){
    cores_needed <- 8 # Limiting #of cores if there are <1000 peaks in total
  }else{
    cores_needed <- total_cores
  }
  doParallel::registerDoParallel(cores = min(total_cores, cores_needed)) # no need for a lot of cores if there aren't that many peaks!
  each_tissue <- floor( min(total_cores, cores_needed) / num_tissues) # Divide cores per tissue and pass onto the foreach loop
  message(paste0("Registering ", min(total_cores, cores_needed), " cores and passing ", each_tissue ," cores per tissue to ", num_tissues ," tissue(s)." ) )
  res_out <- foreach::foreach(tissue = names(qtl_peaks)) %dopar% {
    qtl_mediate(tissue,
                QTL.peaks    = qtl_peaks,
                med_annot    = med_annot,
                QTL.mediator = qtl_mediatior,
                targ_covar   = targ_covar,
                QTL.target   = qtl_target,
                probs        = probs,
                mapDat       = map_dat2,
                cores        = each_tissue
    )
  }
  doParallel::stopImplicitCluster()

  res_list <- list()
  for (i in 1:length(res_out)) {
    # message(i)
    tissue <- res_out[[i]]$tissue
    res_list[[tissue]] <- res_out[[i]]$res_list
  }
  # names(res_list) <- names(qtl_peaks)
  # message(str(res_list))
  if(save %in% c("sr","so")) {
    outfile <- paste0(outdir, "/", med_out)
    saveRDS(res_list, file = outfile)
  }
  if(save %in% c("sr","ro")) {
    return(res_list)
  }
}

qtl_mediate <- function(tissue, QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, mapDat, cores) {
  # res_list <- list()
  # for(tissue in names(QTL.peaks)){
  # message(str(QTL.peaks))
  # message(tissue)
  n.batches <- max(c(round(nrow(QTL.peaks[[tissue]]) / 1000)))
  if( n.batches %in% c(0,1)) n.batches = 2
  nn <- nrow(QTL.peaks[[tissue]])
  ss <- round(seq(0, nn, length.out = n.batches))

  doParallel::registerDoParallel(cores = cores)
  med_res <- foreach::foreach(i = 1:(n.batches - 1)) %dopar% {
    purrr::compact(batchmediate(
      batch = i, QTL.peaks = QTL.peaks[[tissue]],
      med_annot = med_annot[[tissue]],
      QTL.mediator = QTL.mediator[[tissue]],
      targ_covar = targ_covar[[tissue]],
      QTL.target = QTL.target[[tissue]],
      mapDat = mapDat,
      probs = probs[[tissue]],
      ss = ss
    ))
  }
  doParallel::stopImplicitCluster()


  res_list <- do.call("rbind", med_res)
  res_list <- tibble::lst(tissue, res_list)
  return(res_list)
}

batchmediate <- function(batch, z_thres = -2, pos_thres = 10, QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, mapDat, probs, ss, ...) {
  # xx <- parallel::detectCores()
  # yy <- floor(2*xx/(n-1))
  # cl <- parallel::makeCluster(yy)
  # doParallel::registerDoParallel(cl)
  #
  # mutate <- dplyr::mutate
  # select <- dplyr::select
  # filter <- dplyr::filter
  # res_list <- foreach::foreach(batch = 1:(n-1), .combine = "rbind") %dopar% {
  # batch <- n
  med.scan <- list()

  start <- ss[batch] + 1
  end <- ss[batch + 1]
  lod.peaks <- QTL.peaks[start:end, ]
  # cat(sprintf("batch %d: %d-%d\n", batch, start, end))

  for (i in 1:nrow(lod.peaks)) {
    marker <- mapDat |>
      mutate(pos = as.numeric(pos_bp)) |>
      filter(abs(pos - lod.peaks$peak_cM[i]) == min(abs(pos - lod.peaks$peak_cM[i])))
    qtl.chr <- marker$chr
    qtl.pos <- marker$pos_bp / 1e06
    annot <- med_annot |> mutate(middle_point = pos)
    geno <- qtl2::pull_genoprobpos(probs, marker$marker)
    geno <- geno[rownames(geno) %in% rownames(QTL.target), ]
    target <- lod.peaks$phenotype[i]

    med <- intermediate::mediation.scan(
      target = QTL.target[, target, drop = FALSE],
      mediator = QTL.mediator,
      annotation = annot,
      qtl.geno = geno,
      covar = targ_covar,
      verbose = FALSE,
      method = "double-lod-diff"
    )

    med <- med |>
      mutate(
        target_id = lod.peaks$phenotype[i],
        qtl_lod = lod.peaks$lod[i],
        qtl_chr = lod.peaks$peak_chr[i]
      ) |>
      select(
        target_id,
        qtl_lod,
        qtl_chr,
        mediator = symbol,
        mediator_id,
        mediator_chr = chr,
        mediator_midpoint = middle_point,
        LOD
      )

    # return(med)
    med.scan[[i]] <- med
    #   }
  }
  med.scan2 <- dplyr::bind_rows(med.scan)
  return(med.scan2)
  # parallel::stopCluster(cl)
}
