## usethis namespace: start
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import qtl2
## usethis namespace: end
qtl_effects <- function(mapping, peaks, lod_limit, outdir, outfile, n.cores = 4){
  ## Load in data
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
  list2env(tmp_ma,.GlobalEnvp)
  rm(c(tmp_peaks,tmp_map))

  peaks2 <- list()
  for(tissue in names(peaks_list)){
    peaks2[[tissue]] <- peaks_list[[tissue]] %>%
      interp_bp(.) %>%
      dplyr::mutate(interp_bp_peak = ifelse(interp_bp_peak == 3e6, 3000001, interp_bp_peak))
  }

  marker_list <- list()
  for(tissue in names(peaks2)){
    marker_list[[tissue]] <- map_dat2 %>%
      dplyr::select(chrom, pos_bp) %>%
      dplyr::rename(start = pos_bp) %>%
      dplyr::mutate(end = start) %>%
      GenomicRanges::GRanges()
  }

  peaksf <- list()
  for(tissue in names(peaks2)){
    peaksf[[tissue]] <- peaks2[[tissue]] %>%
      dplyr::filter(lod > lod_limit) %>%
      dplyr::arrange(peak_chr, interp_bp_peak)
    query <- peaksf[[tissue]] %>%
      dplyr::select(peak_chr, interp_bp_peak) %>%
      dplyr::rename(chrom = peak_chr, start - interp_bp_peak) %>%
      dplyr::mutate(end = start) %>%
      GenomicRanges::GRanges()
    subject <- marker_list[[tissue]]
    peaksf[[tissue]]$before <- map_dat2$marker[GenomicRanges::follow(query, subject)]
    peaksf[[tissue]]$after <- map_dat2$marker[GenomicRanges::precede(query, subject)]
    peaksf[[tissue]] <- peaksf[[tissue]] %>%
      dplyr::filter(!is.na(before) & !is.na(after))
  }

  effects_blup <- list()
  effects_std <- list()
  foreach::foreach(tissue=names(QTL.peaks)) %dopar% {
    n_peaks <- nrow(peaksf[[tissue]])
    haps <- LETTERS[1:8]
    effects_blupl_temp <- list()
    effects_stdl_temp <- list()
    foreach::foreach(i = 1:n_peaks) %dopar% {
      this_chrom <- peaksf[[tissue]]$peak_chr[i]
      this_markers <- c(peaksf[[tissue]]$before[i], peaksf[[tissue]]$after[i])
      probs_2marker <- subset_probs(qtlprobs[[tissue]], this_chrom, this_markers)
      g <- setNames(list(gmap[[this_chrom]][this_markers]), this_chrom)
      pheno <- peaksf[[tissue]]$phenotype[i]
      out_blup <- qtl2::scan1blup(probs_2marker,
                                  exprZ_list[[tissue]][, pheno, drop = FALSE],
                                  kinship_loco[[tissue]][[this_chrom]],
                                  covar_list[[tissue]],
                                  cores = n.cores)
      out_std <- qtl2::scan1coef(probs_2marker,
                                 exprZ_list[[tissue]][, pheno, drop = FALSE],
                                 kinship_loco[[tissue]][[this_chrom]],
                                 covar_list[[tissue]],
                                 cores = n.cores)
      effects_blupl_temp[[i]] <- colMeans(out_blup[, haps])
      effects_stdl_temp[[i]] <- colMeans(out_std[, haps])
    }
    effects_blup[[tissue]] <- list2DF(effects_blupl_temp[[tissue]])
    effects_std[[tissue]] <- list2DF(effects_stdl_temp[[tissue]])
  }

  peaks <- peaksf

  effects_out <- list(effects_blup, effects_std, peaks)
  saveRDS(effects_out, paste0(outdir,"/",outfile))

  return(effects_out)
}
