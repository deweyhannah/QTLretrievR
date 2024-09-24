#' Prepare and run mediation for a set of QTL peaks
#'
#' @param peaks List of dataframes containing QTL peaks for each tissue
#' @param mapping List of relevant mapping information for overall project
#' @param suggLOD Suggestive LOD to use as filter for mediation. Default is 7.
#' @param outdir String of path to output directory where mediation lists will be saved.
#' @param biomart String pointing to annotations file or annotations object.
#' @param med_out Output file name to save mediation results for later use
#'
#' @return A list containing mediation results for each tissue
#'
#' @export
#'
#' @importFrom dplyr rename filter mutate select
#' @importFrom purrr compact
#' @importFrom BiocParallel bplapply MulticoreParam
#'
run_mediate <- function(peaks, mapping, suggLOD = 7, outdir, biomart, med_out) {
  # mediate_env <- new.env()
  message("load annotations")
  if (is.character(biomart)) {
    biomart <- read.delim(biomart)
  }
  if ("Gene.start..bp." %in% colnames(biomart)) {
    message("renaming biomart columns")
    biomart <- biomart |>
      dplyr::rename(gene = Gene.stable.ID, symbol = MGI.symbol, start = Gene.start..bp., end = Gene.end..bp., chr = Chromosome.scaffold.name)
  }
  if ("gene.id" %in% colnames(biomart)) {
    colnames(biomart)[which(colnames(biomart) == "gene.id")] <- "gene"
  }
  if ("end" %in% colnames(biomart)) {
    colnames(biomart)[which(colnames(biomart) == "end")] <- "stop"
  }
  if ("chr" %in% colnames(biomart)) {
    colnames(biomart)[which(colnames(biomart) == "chr")] <- "chromosome"
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

  biomart_list <- list()
  for (tissue in names(peaks_list)) {
    biomart_list[[tissue]] <- biomart |>
      dplyr::filter(gene %in% rownames(expr_list[[tissue]])) |>
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
      dplyr::filter(target_id %in% biomart_list[[tissue]]$gene)
    qtl_target[[tissue]] <- exprZ_list[[tissue]][, biomart_list[[tissue]]$gene]
  }

  message("filtered peaks")

  ## Run annotations
  targ_annot <- list()
  for (tissue in names(peaks_list)) {
    # message(paste0(names(biomart_list[[tissue]]), sep = " "))
    targ_annot[[tissue]] <- biomart_list[[tissue]] |>
      dplyr::mutate(target_id = gene, chrom = chromosome) |>
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
    qtl_mediatior[[tissue]] <- exprZ_list[[tissue]][, biomart_list[[tissue]]$gene]
    qtl_mediatior[[tissue]] <- qtl_mediatior[[tissue]][, med_annot[[tissue]]$gene, drop = FALSE]
    qtl_target[[tissue]] <- qtl_target[[tissue]][, targ_annot[[tissue]]$target_id, drop = FALSE]
    qtl_peaks[[tissue]] <- qtl_peaks[[tissue]] |>
      dplyr::filter(target_id %in% targ_annot[[tissue]]$target_id)
  }
  med_covar <- covar_list

  message("running mediation")

  each_tissue <- floor( as.numeric(parallelly::availableCores()) / length(names(exprZ_list)))
  doParallel::registerDoParallel(cores = each_tissue)
  res_out <- foreach::foreach(tissue = names(qtl_peaks)) %dopar% {
    qtl_mediate(tissue,
                QTL.peaks = qtl_peaks, med_annot = med_annot, QTL.mediator = qtl_mediatior,
                targ_covar = targ_covar, QTL.target = qtl_target, probs = probs,
                mapDat = map_dat2
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
  outfile <- paste0(outdir, "/", med_out)
  saveRDS(res_list, file = outfile)
  return(res_list)
}

qtl_mediate <- function(tissue, QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, mapDat) {
  # res_list <- list()
  # for(tissue in names(QTL.peaks)){
  # message(str(QTL.peaks))
  # message(tissue)
  n.batches <- max(c(round(nrow(QTL.peaks[[tissue]]) / 1000)))
  nn <- nrow(QTL.peaks[[tissue]])
  ss <- round(seq(0, nn, length.out = n.batches))

  # med.scans <- foreach::foreach(i = 1:(n.batches - 1), .combine = "rbind") %dopar% {
  med_res <- BiocParallel::bplapply(1:(n.batches - 1), function(i) {
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
  })
  # }
  # outfile <- paste0(outdir, "/", med_out)
  # saveRDS(res_list, file = outfile)

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

    return(med)

    #   }
  }
  # parallel::stopCluster(cl)
}
