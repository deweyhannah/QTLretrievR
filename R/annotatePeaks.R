#' Interpolate peak positions and attach annotations
#'
#' @param map map_list from [mapQTL].
#' @param peaks Un-annotated peaks.
#' @param annots Annotations file. Path to location or object.
#' @param localRange is defined as "local". Default is 10e6.
#'
#' @return List of annotated peaks for each tissue.
#'
#' @export
#'
#' @importFrom utils read.delim
#' @importFrom dplyr left_join mutate select rename
#'
#'
annotatePeaks <- function(map, peaks, annots, localRange = 10e6) {
  ## Get annots columnnames to wanted format from base download
  if (is.character(annots)) {
    annots <- read.delim(annots)
  }
  if ("Gene.start..bp." %in% colnames(annots)) {
    message("renaming annotation columns")
    annots <- annots |>
      dplyr::rename(id = Gene.stable.ID, symbol = MGI.symbol, start = Gene.start..bp., end = Gene.end..bp., chr = Chromosome.scaffold.name)
  }
  if ("gene.id" %in% colnames(annots)) {
    colnames(annots)[which(colnames(annots) == "gene.id")] <- "id"
  }

  ## Add midpoint of gene to annotations and pare down the columns to those we need
  annots <- annots |>
    dplyr::mutate(midpoint = (start + end) / 2) |>
    dplyr::select(id, symbol, chr, start, end, midpoint)

  ## Interpolate the physical bp locations from centimorgans join to annotationss
  peaks_pmap <- list()
  for (tissue in names(peaks)) {
    # message(paste0(colnames(peaks[[tissue]]), sep = " "))
    peaks_pmap[[tissue]] <- peaks[[tissue]] |>
      # interp_bp(genmap = map$gmap, physmap = map$pmap) |>
      # dplyr::mutate(phenotype = gsub("_.*", "", phenotype)) |>
      dplyr::left_join(annots, by = c("phenotype" = "id"))
  }

  ## Annotate for locality
  annotated_peaks <- list()
  for (tissue in names(peaks_pmap)) {
    annotated_peaks[[tissue]] <- peaks_pmap[[tissue]] |>
      dplyr::mutate(local = ifelse(abs(midpoint - peak_bp) < localRange & chr == peak_chr, 1, 0))
  }

  return(annotated_peaks)
}
