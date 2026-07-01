#' Prepare and run mediation for delta G x E peaks
#'
#' @param peaks List of dataframes containing delta G x E peaks (name 'delta')
#' @param mapping Mapping list from `gxeQTL`. Needs to contain mapping information
#'  for at minimum the delta condition, and the condition being used for mediation
#' @param by String indicating the name of the condition to be used for mediation.
#'  Default 'delta'.
#' @param sigLOD Significant LOD threshold to use to filter phenotypes for
#'  mediation. Default is 7.5
#' @param annots Annotations file. Contains mapping information for phenotypes.
#'  Dataframe, or tsv. Columns must include "id", "symbol", "start", "end".
#' @param outdir Directory to save output files. Default is `NULL`.
#' @param med_out String indicating the name of the output file containing
#'  mediation results for mediation within a phenotype. This file will be
#'   saved in `.rds` format and used for downstream analysis and visualization.
#'    Should end in `.rds`. Default is "`delta_mediation.rds`"
#' @param total_cores Number of available cores to use for parallelization.
#'  Default is `NULL`.
#' @param save Indicates object return/save behavior. One of
#'  `c("sr", "so", "ro")`; save & return, save only, return only.
#'   Default is "sr".
#' @param distOnly Logical. Mediate only the distal peaks? Default is `TRUE`.
#' @param hsOnly Logical. Mediate only on peaks identified within a hotspot.
#'  Default is `FALSE`.
#' @param BPPARAM BiocParallel Parameter
#'
#' @return A list containing mediation results for the peaks called from the delta G x E method.
#'  Mediated by 'delta', '<ctrl>', or '<env>'.
#'
#' @export
#'
#' @importFrom dplyr rename filter mutate select
#' @importFrom purrr compact
#' @importFrom BiocParallel SerialParam MulticoreParam bpnworkers
#'
gxe_modiFinder <- function(peaks, mapping, by="delta", sigLOD = 7.5, annots, outdir = NULL,
                       med_out = "gxe_mediation.rds", total_cores = NULL,
                       save = "sr", distOnly = TRUE, hsOnly = FALSE,
                       BPPARAM = BiocParallel::SerialParam()) {
  ## Check save conflicts
  if (save %in% c("sr", "so")) {
    if (is.null(outdir)) {
      stop("Requested Save. No output directory provided, and no default.")
    }
  }

  message(paste0("Mediating G x E peaks by ", by, " rankZ transformed expression. If this was not intended please check the 'by' argument."))

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
    # tmp_map <- check_data(mapping)
  }
  if (is.character(peaks)) {
    peaks_list <- check_data(peaks, type = "peaks")
    # tmp_map <- check_data(mapping)
  }

  if (!is.data.frame(peaks_list[[1]])) {
    stop("The elements of 'peaks_list' must be data frames.")
  }

  message("data checked")

  ## clean up
  if (!is.null(mapping)) {
    exprZ_gxe <- mapping$exprZ_list
    covar_list <- mapping$covar_list
    expr_list <- mapping$expr_list
    gmap <- mapping$gmap
    kinship_loco <- mapping$kinship_loco
    map_dat2 <- mapping$map_dat2
    pmap <- mapping$pmap
    qtlprobs <- mapping$qtlprobs
    tissue_samp <- mapping$tissue_samp
  }
  rm(mapping)

  annots_list <- list()
  for (tissue in names(peaks_list)) {
    annots_list[[tissue]] <- annots |>
      dplyr::filter(id %in% colnames(exprZ_gxe[[tissue]])) |>
      dplyr::mutate(midpoint = (start + stop) / 2)
  }

  if (hsOnly) {
    tbands <- transbands(map_dat = map_dat2,
                         peaks   = peaks_list,
                         sigLOD  = sigLOD,
                         psave   = FALSE)

    peaks_og <- peaks_list
    for (tissue in names(peaks_list)) {
      peaks_list[[tissue]] <- filter_peaks(peaks_list[[tissue]],
                                           tbands$bands.rna[[tissue]])
    }
  }

  ## Identify tissues and targets
  qtl_peaks <- list()
  qtl_target <- list()
  for (tissue in names(peaks_list)) {
    qtl_peaks[[tissue]] <- peaks_list[[tissue]] |>
      dplyr::filter(lod > sigLOD) |>
      # interp_bp(.) |>
      dplyr::mutate(phenotype = gsub("_.*", "", phenotype)) |>
      dplyr::mutate(target_id = phenotype) |>
      dplyr::filter(target_id %in% annots_list[[tissue]]$id)
    qtl_target[[tissue]] <- exprZ_gxe[[tissue]][, annots_list[[tissue]]$id]
  }

  message("filtered peaks")

  ## Run annotations
  targ_annot <- list()
  for (tissue in names(peaks_list)) {
    targ_annot[[tissue]] <- annots_list[[tissue]] |>
      dplyr::mutate(target_id = id, chrom = chromosome) |>
      dplyr::mutate(chrom = ifelse(chromosome == "MT", "M", chromosome)) |>
      dplyr::filter(!is.na(chrom)) |>
      dplyr::mutate(chr = chrom, pos = abs(stop + start) / 2)
  }
  targ_covar <- covar_list
  probs <- qtlprobs
  kinship <- kinship_loco
  # return(targ_covar)

  ## mediator's annotation
  med_annot <- annots |>
    dplyr::filter(id %in% colnames(exprZ_gxe[[by]])) |>
    dplyr::mutate(mediator_id = id, chrom = chromosome) |>
    dplyr::mutate(chrom = ifelse(chromosome == "MT", "M", chromosome)) |>
    dplyr::filter(!is.na(chrom)) |>
    dplyr::mutate(chr = chrom, pos = abs(stop + start) / 2)

  ## mediator expression
  qtl_mediator <- list()
  for (tissue in names(peaks_list)) {
    qtl_mediator[[tissue]] <- exprZ_gxe[[by]][, med_annot$mediator_id]
  }

  ## mediatior covariates
  med_covar <- covar_list


  message("running mediation")

  # available_cores <- get_cores()
  # if( is.null(total_cores)) total_cores <- available_cores
  # if( total_cores > available_cores) total_cores <- available_cores
  # get the maximum number of peaks
  max_peaks <- max(vapply(qtl_peaks, nrow, integer(1)))
  num_tissues <-  length(names(qtl_peaks)) # number of tissues
  if( max_peaks < 1000){
    cores_needed <- 8 # Limiting # of cores if there are <1000 peaks in total
  }else{
    cores_needed <- workers
  }

  res_out <- lapply(names(qtl_peaks), function(tissue) {
    multi_qtl_mediate(
      tissue        = tissue,
      QTL.peaks     = qtl_peaks,
      med_annot     = med_annot,
      QTL.mediator  = qtl_mediator,
      targ_covar    = targ_covar,
      QTL.target    = qtl_target,
      probs         = probs,
      mapDat        = map_dat2,
      pmap          = pmap,
      BPPARAM       = BPPARAM
    )})

  res_list <- list()
  for (i in seq_len(length(res_out))) {
    tissue <- res_out[[i]]$tissue
    res_list[[tissue]] <- res_out[[i]]$res_list
    res_list[[tissue]] <- res_list[[tissue]] |>
      dplyr::left_join(annots |> dplyr::select(target_id = id, target = symbol))
  }
  if(save %in% c("sr","so")) {
    outfile <- paste0(outdir, "/", med_out)
    saveRDS(res_list, file = outfile)
  }
  if(save %in% c("sr","ro")) {
    return(res_list)
  }
}
