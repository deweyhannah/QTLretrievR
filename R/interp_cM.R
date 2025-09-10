#' interp_cM
#' @description
#' Interpolate bp peak positions to cM based on marker maps.
#'
#'
#' @param df A peak dataframe with the peak position in bp.
#' @param genmap Genetic map - found in `maps_list` (gmap).
#' @param physmap Physical map - found in `maps_list` (pmap).
#'
#' @return Peak dataframe with interpolated bp column appended
#'
#' @export
#'
interp_cM <- function(df, genmap, physmap) {
  if(peak_bp %in% colnames(df)) {
    message("genomic peak locations already calculated, did you mean to
            use interp_bp?")
    return(df)
  }
  chroms <- c(as.character(1:19), "X")
  df <- dplyr::arrange(df, peak_chr, peak_bp)
  peak_ppos <- select(df, peak_chr, peak_bp)
  chr <- peak_ppos$peak_chr
  f <- factor(chr, chroms)
  peak_pcoord_list <- split(peak_ppos$peak_bp, f)
  peak_gcoord_list <- qtl2::interp_map(peak_pcoord_list, physmap, genmap)
  df$peak_bp <- unsplit(peak_gcoord_list, f)
  df
}
