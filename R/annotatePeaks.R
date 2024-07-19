## usethis namespace: start
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
## usethis namespace: end
annotatePeaks <- function(map, peaks, biomart, localRange = 10e6){
  ## Get biomart columnnames to wanted format from base download
  if("Gene.start..bp." %in% colnames(biomart)){
    message("renaming biomart columns")
    biomart <- biomart %>%
      dplyr::rename(gene = Gene.stable.ID, symbol = MGI.symbol, start = Gene.start..bp., end = Gene.end..bp., chr = Chromosome.scaffold.name)
  }

  ## Add midpoint of gene to biomart and pare down the columns to those we need
  biomart <- biomart %>%
    dplyr::mutate(midpoint = (start+end)/2) %>%
    dplyr::select(gene, symbol, chr, start, end, midpoint)

  ## Interpolate the physical bp locations from centimorgans join to biomart
  peaks_pmap <- list()
  for(tissue in names(peaks)){
    peaks_pmap[[tissue]] <- peaks[[tissue]] %>%
      interp_bp(.) %>%
      dplyr::mutate(phenotype = gsub('_.*','', phenotype)) %>%
      dplyr::left_join(biomart, by = join_by(phenotype == gene))
  }

  ## Annotate for locality
  annotated_peaks <- list()
  for(tissue in names(peaks_pmap)){
    annotated_peaks[[tissue]] <- peaks_pmap[[tissue]] %>%
      dplyr::mutate(local = ifelse(abs(midpoint - interp_bp_peak) < localRange & chr == peak_chr, 1, 0))
  }

  return(annotated_peaks)
}
