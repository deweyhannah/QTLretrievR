## usethis namespace: start
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom purrr compact
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
## usethis namespace: end
prep_mediate <- function(peaks, mapping, peak_lim = 7, outdir, biomart, effects_out){
  if("Gene.start..bp." %in% colnames(biomart)){
    message("renaming biomart columns")
    biomart <- biomart %>%
      dplyr::rename(gene = Gene.stable.ID, symbol = MGI.symbol, start = Gene.start..bp., end = Gene.end..bp., chr = Chromosome.scaffold.name)
  }
  if((is.character(peaks) & is.list(mapping)) | (is.list(peaks) & is.character(mapping))){
    stop("Peaks and mapping must both direct to an RDS file or be lists")
  }

  ## Load and organize relevant data
  if(is.list(peaks)){
    tmp_peaks <- check_data(peaks)
    tmp_map <- check_data(mapping)
  }
  if(is.character(peaks)){
    tmp_peaks <- check_data(paste0(outdir,"/",peaks))
    tmp_map <- check_data(paste0(outdir,"/",mapping))
  }
  list2env(tmp_peaks,.GlobalEnv)
  list2env(tmp_map,.GlobalEnv)
  rm(c(tmp_peaks,tmp_map))

  biomart_list <- list()
  for(tissue in names(peaks_list)){
    biomart_list[[tissue]] <- biomart %>%
      dplyr::filter(gene %in% rownames(expr_list[[tissue]])) %>%
      dplyr::mutate(midpoint = (start+stop)/2)
  }

  ## Identify tissues and targets
  qtl_peaks <- list()
  qtl_target <- list()
  for(tissue in names(peaks_list)){
    qtl_peaks[[tissue]] <- peaks_list[[tissue]] %>%
      dplyr::filter(lod > peak_lim) %>%
      interp_bp(.) %>%
      dplyr::mutate(phenotype = gsub('_.*',"",phenotype)) %>%
      dplyr::mutate(target_id = phenotype) %>%
      dplyr::filter(target_id %in% biomart_list[[tissue]]$gene)
    qtl_target[[tissue]] <- exprZ_list[[tissue]][,biomart_list[[tissue]]$gene]
  }

  ## Run annotations
  targ_annot <- list()
  for(tissue in names(peaks_list)){
    targ_annot[[tissue]] <- biomart_list[[tissue]] %>%
      dplyr::mutate(target_id = gene, chrom = chromosome) %>%
      dplyr::mutate(chrom = ifelse(chromosome == "MT", "M", chromosome)) %>%
      dplyr::filter(!is.na(chrom)) %>%
      dplyr::mutate(chr = chrom, pos = abs(stop+start)/2)
  }
  targ_covar <- covar_list
  probs <- qtlprobs
  kinship <- kinship_loco

  qtl_mediatior <- list()
  med_annot <- list()
  for(tissue in names(exprZlist)){
    med_annot[[tissue]] <- targ_annot[[tissue]] %>%
      dplyr::rename(mediator_id = target_id)
    qtl_mediatior[[tissue]] <- exprZ_list[[tissue]][,biomart_list[[tissue]]$gene]
    qtl_mediatior[[tissue]] <- qtl_mediatior[[tissue]][,med_annot[[tissue]]$gene, drop = FALSE]
    qtl_target[[tissue]] <- qtl_target[[tissue]][,targ_annot[[tissue]]$target_id, drop = FALSE]
    qtl_peaks[[tissue]] <- qtl_peaks[[tissue]] %>%
      dplyr::filter(target_id %in% targ_annot[[tissue]]$target_id)
  }
  med_covar <- covar_list

  res_list <- qtl_mediate(QTL.peaks = qtl_peaks, med_annot = med_annot, QTL.mediator = qtl_mediatior,
                          targ_covar = targ_covar, QTL.target = qtl_target, probs = probs,
                          effects_out = effects_out, outdir = outdir)
  return(res_list)
}

qtl_mediate <- function(QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, effects_out, outdir){
  med.scans <- list()
  results <- list()
  res_list <- list()
  foreach::foreach(tissue=names(QTL.peaks)) %dopar% {
    n.batches <- max(c(round(nrow(QTL.peaks)/1000)))
    nn <- nrow(QTL.peaks)
    ss <- round(seq(0, nn, length.out = n.batches))
    message(paste0("Mediating ", nn, " eQTL peaks. Running in ", n.batches - 1, " batches"))
    med.scans[[tissue]] <- list()
    results[[tissue]] <- list()
    foreach::foreach(i = 1:(n.batches - 1)) %dopar% {
      med.scans[[tissue]][[i]] <- purrr::compact(batchmediate(n = i, QTL.peaks = QTL.peaks[[tissue]],
                                                              med_annot = med_annot[[tissue]],
                                                              QTL.mediator = QTL.mediator[[tissue]],
                                                              targ_covar = targ_covar[[tissue]],
                                                              QTL.target = QTL.target[[tissue]],
                                                              probs = probs[[tissue]]))
      results[[tissue]][[i]] <- list2DF(med.scans[[tissue]][[i]])
    }
    res_list[[tissue]] <- list2DF(results[[tissue]])
  }
  outfile <- paste0(outdir, "/", med_out)
  saveRDS(res_list, file = outfile)

  return(res_list)
}
