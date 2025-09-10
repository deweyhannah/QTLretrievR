#' Prepare and run between phenotype mediation for a set of QTL peaks with
#'  provided expression data.
#'
#' @description
#' Mediation scan between two different phenotypes. Provide peaks and mapping
#' for the phenotype to mediate, and the rank Z-transformed phenotype
#' quantification to use as mediators.
#'
#'
#' @param peaks List of dataframes containing QTL peaks for each tissue, or
#'  full path pointing to saved peaks object.
#' @param mapping Mapping list from `mapQTL`, or full path to `.rds`
#'  containing one.
#' @param exprZ rankZ transformed expression to use for the mediation.
#' @param suggLOD Significant LOD threshold to use to filter phenotypes for
#'  mediation. Default is 7.5
#' @param annots Annotations file. Contains mapping information for phenotypes.
#'  Dataframe, or tsv. Columns must include "id", "symbol", "start", "end".
#' @param outdir Directory to save output files. Default is `NULL`.
#' @param med_out String indicating the name of the output file containing
#' mediation results for mediation between phenotypes. This file will be saved
#'  in `.rds` format and used for downstream analysis and visualization. Should
#'   end in `.rds`. Default is "`multi_pheno_mediation.rds`"
#' @param total_cores Number of available cores to use for parallelization.
#'  Default is `NULL`.
#' @param save Indicates object return/save behavior. One of
#'  `c("sr", "so", "ro")`; save & return, save only, return only.
#'   Default is "sr".
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
multi_modiFinder <- function(peaks, mapping, exprZ, suggLOD = 7, annots,
                             outdir = NULL,
                             med_out = "multi_pheno_mediation.rds",
                             total_cores = NULL, save = "sr") {

  ## Check save conflicts
  if (save %in% c("sr", "so")) {
    if (is.null(outdir)) {
      stop("Requested Save. No output directory provided, and no default.")
    }
  }

  message("load annotations")
  if (is.character(annots)) {
    annots <- read.delim(annots)
  }
  if ("Gene.start..bp." %in% colnames(annots)) {
    message("renaming annots columns")
    annots <- annots |>
      dplyr::rename(id = Gene.stable.ID, symbol = MGI.symbol,
                    start = Gene.start..bp., end = Gene.end..bp.,
                    chr = Chromosome.scaffold.name)
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
  if ((is.character(peaks) & is.list(mapping)) | (is.list(peaks) &
                                                  is.character(mapping))) {
    stop("Peaks and mapping must both direct to an RDS file or be lists")
  }

  ## Load and organize relevant data
  if (is.list(peaks)) {
    peaks_list <- check_data(peaks, type = "peaks")
    tmp_map <- check_data(mapping)
  }
  if (is.character(peaks)) {
    peaks_list <- check_data(peaks, type = "peaks")
    tmp_map <- check_data(mapping)
  }

  if (!is.data.frame(peaks_list[[1]])) {
    stop("The elements of 'peaks_list' must be data frames.")
  }

  message("data checked")

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

  common_samples <- list()
  for (tissue in names(peaks_list)) {
    common_samples[[tissue]] <- intersect(rownames(exprZ),
                                          rownames(exprZ_list[[tissue]]))
    message(paste0("Number of common samples for ", tissue, ": ",
                   length(common_samples[[tissue]])))
  }

  ## Identify tissues and targets
  qtl_peaks <- list()
  qtl_target <- list()
  for (tissue in names(peaks_list)) {
    qtl_peaks[[tissue]] <- peaks_list[[tissue]] |>
      dplyr::filter(lod > suggLOD) |>
      dplyr::mutate(target_id = phenotype)
    qtl_target[[tissue]] <- exprZ_list[[tissue]][common_samples[[tissue]],,
                                                 drop = FALSE]
  }

  message("filtered peaks")

  targ_covar <- list()
  for(tissue in names(qtl_target)) {
    message(tissue)
    targ_covar[[tissue]] <- covar_list[[tissue]][common_samples[[tissue]], ,
                                                 drop = FALSE]
  }
  probs <- qtlprobs
  kinship <- kinship_loco


  ## mediator's annotation
  med_annot <- annots |>
    dplyr::filter(id %in% colnames(exprZ)) |>
    dplyr::mutate(mediator_id = id, chrom = chromosome) |>
    dplyr::mutate(chrom = ifelse(chromosome == "MT", "M", chromosome)) |>
    dplyr::filter(!is.na(chrom)) |>
    dplyr::mutate(chr = chrom, pos = abs(stop + start) / 2)

  ## mediator expression
  qtl_mediator <- list()
  for (tissue in names(peaks_list)) {
    qtl_mediator[[tissue]] <- exprZ[common_samples[[tissue]],
                                    med_annot$mediator_id]
  }

  ## mediatior covariates
  med_covar <- covar_list

  message("running mediation")

  available_cores <- get_cores()
  if( is.null(total_cores)) total_cores <- available_cores
  if( total_cores > available_cores) total_cores <- available_cores
  ## get the maximum number of peaks
  max_peaks <- max(vapply(qtl_peaks, nrow, integer(1)))
  num_tissues <-  length(names(qtl_peaks)) # number of tissues
  if( max_peaks < 1000){
    cores_needed <- 8 ## Limiting # of cores if there are <1000 peaks in total
  }else{
    cores_needed <- total_cores
  }

  doParallel::registerDoParallel(cores = min(total_cores, cores_needed))
  each_tissue <- floor( min(total_cores, cores_needed) / num_tissues)
  message(paste0("Registering ", min(total_cores, cores_needed),
                 " cores and passing ", each_tissue ," cores per tissue to ",
                 num_tissues ," tissue(s)." ) )
  res_out <- foreach::foreach(tissue = names(qtl_peaks)) %dopar% {
    multi_qtl_mediate(tissue,
                QTL.peaks    = qtl_peaks,
                med_annot    = med_annot,
                QTL.mediator = qtl_mediator,
                targ_covar   = targ_covar,
                QTL.target   = qtl_target,
                probs        = probs,
                mapDat       = map_dat2,
                cores        = each_tissue,
                pmap         = pmap
    )
  }
  doParallel::stopImplicitCluster()

  res_list <- list()
  for (i in seq_len(length(res_out))) {
    tissue <- res_out[[i]]$tissue
    res_list[[tissue]] <- res_out[[i]]$res_list
  }

  if(save %in% c("sr","so")) {
    outfile <- paste0(outdir, "/", med_out)
    saveRDS(res_list, file = outfile)
  }

  ## Return both mapping and peak results as a named list
  if(save %in% c("sr","ro")) {
    return(res_list)
  }
}

multi_qtl_mediate <- function(tissue, QTL.peaks, med_annot, QTL.mediator,
                              targ_covar, QTL.target, probs, mapDat, cores,
                              pmap) {

  n.batches <- max(c(round(nrow(QTL.peaks[[tissue]]) / 1000)))
  if( n.batches %in% c(0,1)) n.batches <- 2
  nn <- nrow(QTL.peaks[[tissue]])
  ss <- round(seq(0, nn, length.out = n.batches))

  doParallel::registerDoParallel(cores = cores)
  med_res <- foreach::foreach(i = seq_len(n.batches - 1)) %dopar% {
    purrr::compact(batchmediate(
      batch        = i,
      QTL.peaks    = QTL.peaks[[tissue]],
      med_annot    = med_annot,
      QTL.mediator = QTL.mediator[[tissue]],
      targ_covar   = targ_covar[[tissue]],
      QTL.target   = QTL.target[[tissue]],
      mapDat       = mapDat,
      probs        = probs[[tissue]],
      ss           = ss,
      pmap         = pmap
    ))
  }
  doParallel::stopImplicitCluster()


  res_list <- do.call("rbind", med_res)
  res_list <- tibble::lst(tissue, res_list)
  return(res_list)
}
