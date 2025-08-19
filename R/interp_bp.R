#' interp_bp
#' @description
#' Interpolate cM peak positions to bp based on marker maps
#'
#'
#' @param df A peak dataframe with the peak position in cM
#' @param genmap Genetic map - found in `maps_list` (gmap)
#' @param physmap Physical map - found in `maps_list` (pmap)
#'
#' @return Peak dataframe with interpolated bp column appended
#'
#' @export
#'
interp_bp <- function(df, genmap, physmap) {
  if("peak_bp" %in% colnames(df)) {
    message("physical peak locations already calculated, did you mean to use interp_cM?")
    return(df)
  }
  chroms <- c(as.character(1:19), "X")
  df <- dplyr::arrange(df, peak_chr, peak_cM)
  peak_gpos <- select(df, peak_chr, peak_cM)
  chr <- peak_gpos$peak_chr
  f <- factor(chr, chroms)
  peak_gcoord_list <- split(peak_gpos$peak_cM, f)
  peak_pcoord_list <- qtl2::interp_map(peak_gcoord_list, genmap, physmap)
  df$peak_bp <- unsplit(peak_pcoord_list, f)
  df
}
