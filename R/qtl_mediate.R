## usethis namespace: start
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom purrr compact
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
## usethis namespace: end
prep_mediate <- function(peaks, mapping, peak_lim = 7, outdir, biomart, effects_out, n.cores = 4){
  # mediate_env <- new.env()
  message("load annotations")
  if(is.character(biomart)){
    biomart <- read.delim(biomart)
  }
  if("Gene.start..bp." %in% colnames(biomart)){
    message("renaming biomart columns")
    biomart <- biomart %>%
      dplyr::rename(gene = Gene.stable.ID, symbol = MGI.symbol, start = Gene.start..bp., end = Gene.end..bp., chr = Chromosome.scaffold.name)
  }
  if("gene.id" %in% colnames(biomart)){
    colnames(biomart)[which(colnames(biomart)=="gene.id")] <- "gene"
  }
  if("end" %in% colnames(biomart)){
    colnames(biomart)[which(colnames(biomart)=="end")] <- "stop"
  }
  if("chr" %in% colnames(biomart)){
    colnames(biomart)[which(colnames(biomart)=="chr")] <- "chromosome"
  }

  message("checking peaks and mapping")
  if((is.character(peaks) & is.list(mapping)) | (is.list(peaks) & is.character(mapping))){
    stop("Peaks and mapping must both direct to an RDS file or be lists")
  }

  ## Load and organize relevant data
  if(is.list(peaks)){
    peaks_list <- check_data(peaks, type = "peaks")
    tmp_map <- check_data(mapping)
    # if(!is.null(tmp_map)){
    #   list2env(mapping, .GlobalEnv)
    # }
  }
  if(is.character(peaks)){
    peaks_list <- check_data(paste0(outdir,"/",peaks), type = "peaks")
    tmp_map <- check_data(paste0(outdir,"/",mapping))
    # message(paste0(names(tmp_map), sep = " "))
  }
  message("data checked")

  ## clean up
  # list2env(tmp_peaks,.GlobalEnv)
  # list2env(tmp_map,.GlobalEnv)
  # rm(c(tmp_peaks,tmp_map))
  if(!is.null(tmp_map)){
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
  for(tissue in names(peaks_list)){
    biomart_list[[tissue]] <- biomart %>%
      dplyr::filter(gene %in% rownames(expr_list[[tissue]])) %>%
      dplyr::mutate(midpoint = (start+stop)/2)
  }

  # message(paste0(ls(), sep = " "))

  ## Identify tissues and targets
  qtl_peaks <- list()
  qtl_target <- list()
  for(tissue in names(peaks_list)){
    qtl_peaks[[tissue]] <- peaks_list[[tissue]] %>%
      dplyr::filter(lod > peak_lim) %>%
      # interp_bp(.) %>%
      dplyr::mutate(phenotype = gsub('_.*',"",phenotype)) %>%
      dplyr::mutate(target_id = phenotype) %>%
      dplyr::filter(target_id %in% biomart_list[[tissue]]$gene)
    qtl_target[[tissue]] <- exprZ_list[[tissue]][,biomart_list[[tissue]]$gene]
  }

  message("filtered peaks")

  ## Run annotations
  targ_annot <- list()
  for(tissue in names(peaks_list)){
    # message(paste0(names(biomart_list[[tissue]]), sep = " "))
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
  for(tissue in names(exprZ_list)){
    med_annot[[tissue]] <- targ_annot[[tissue]] %>%
      dplyr::rename(mediator_id = target_id)
    qtl_mediatior[[tissue]] <- exprZ_list[[tissue]][,biomart_list[[tissue]]$gene]
    qtl_mediatior[[tissue]] <- qtl_mediatior[[tissue]][,med_annot[[tissue]]$gene, drop = FALSE]
    qtl_target[[tissue]] <- qtl_target[[tissue]][,targ_annot[[tissue]]$target_id, drop = FALSE]
    qtl_peaks[[tissue]] <- qtl_peaks[[tissue]] %>%
      dplyr::filter(target_id %in% targ_annot[[tissue]]$target_id)
  }
  med_covar <- covar_list

  message("running mediation")
  res_list <- qtl_mediate(QTL.peaks = qtl_peaks, med_annot = med_annot, QTL.mediator = qtl_mediatior,
                          targ_covar = targ_covar, QTL.target = qtl_target, probs = probs,
                          effects_out = effects_out, outdir = outdir, n.cores = n.cores)
  return(res_list)
}

qtl_mediate <- function(QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, effects_out, outdir, n.cores){
  xx <- parallel::detectCores()
  yy <- floor(xx/length(names(QTL.peaks)))
  cl <- parallel::makeCluster(yy)
  doParallel::registerDoParallel(cl)
  for(tissue in names(QTL.peaks)){
    n.batches <- max(c(round(nrow(QTL.peaks[[tissue]])/1000)))
    nn <- nrow(QTL.peaks[[tissue]])
    message(paste0("Mediating ", nn, " eQTL peaks. Running in ", n.batches - 1, " batches"))
    ss <- round(seq(0, nn, length.out = n.batches))
    xx <- parallel::detectCores()
    yy <- floor(xx/n.cores)
    cl <- parallel::makeCluster(yy)
    doParallel::registerDoParallel(cl)
    # message(paste0("Mediating ", nn, " eQTL peaks. Running in ", n.batches - 1, " batches"))
    # cl2 <- parallel::makeCluster(floor(yy/n.cores))
    doParallel::registerDoParallel(cl)
    med.scans <- foreach::foreach(i = 1:(n.batches - 1), .combine = "rbind") %dopar% {
      scans <- batchmediate(n = i, QTL.peaks = QTL.peaks[[tissue]],
                   med_annot = med_annot[[tissue]],
                   QTL.mediator = QTL.mediator[[tissue]],
                   targ_covar = targ_covar[[tissue]],
                   QTL.target = QTL.target[[tissue]],
                   probs = probs[[tissue]])
      scans %>% purrr::compact()
    }
    parallel::stopCluster(cl)
    res_list[[tissue]] <- med.scans
    # res_list[[tissue]] <- list2DF(results[[tissue]])
  }
  # parallel::stopCluster(cl)
  # names(res_list) <- names(QTL.peaks)
  outfile <- paste0(outdir, "/", med_out)
  saveRDS(res_list, file = outfile)

  return(res_list)
}

batchmediate <- function( n, z_thres = -2,  pos_thres = 10, QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, ...){

  med.scan <- list()

  start <- ss[n]+1
  end   <- ss[n+1]
  lod.peaks <- QTL.peaks[start:end,]
  cat(sprintf("batch %d: %d-%d\n", n, start, end))

  for(i in 1:nrow(lod.peaks) ){

    marker    <- map_dat2 %>%
      dplyr::mutate(pos=as.numeric(pos_bp)) %>%
      dplyr::filter( abs(pos - lod.peaks$peak_cM[i]) == min(abs(pos - lod.peaks$peak_cM[i])))
    qtl.chr   <- marker$chr
    qtl.pos   <- marker$pos_bp/1e06
    annot     <- med_annot %>% mutate( middle_point = pos)
    geno      <- qtl2::pull_genoprobpos(probs,marker$marker)
    geno      <- geno[rownames(geno) %in% rownames(QTL.target),]
    target    <- lod.peaks$phenotype[i]

    med <- intermediate::mediation.scan(target     = QTL.target[, target, drop = FALSE],
                                        mediator   = QTL.mediator,
                                        annotation = annot,
                                        qtl.geno   = geno,
                                        covar      = targ_covar,
                                        verbose    = FALSE,
                                        method     = "double-lod-diff")

    med <- med %>%
      mutate(
        target_id = lod.peaks$phenotype[i],
        qtl_lod = lod.peaks$lod[i],
        qtl_chr = lod.peaks$peak_chr[i]
      ) %>%
      select(
        target_id,
        qtl_lod,
        qtl_chr,
        mediator = symbol,
        mediator_id ,
        mediator_chr = chr,
        mediator_midpoint = middle_point,
        LOD)

    med.scan[[i]] <- med

    print(i)
  }
}
