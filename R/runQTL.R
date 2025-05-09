#' Wrapper function to generate mapping, peaks, mediation, and effects data
#'
#' @param geno_out Output file name to save gbrs interpolated genoprobs. Default is "gbrs_interpolated_genoprobs.rds".
#' @param peaks_out Output file name to save peaks. Default is "mm39_peaks.rds".
#' @param map_out Output file name to save mapping. Default is "mm39_mapping.rds".
#' @param med_out Output file name to save mediation. Default is "mm39_mediation_res.rds".
#' @param effects_out Output file name to save effects. Default is "mm39_effects.rds".
#' @param outdir Output directory to save results.
#' @param gbrs_fileLoc Path to GBRS interpolated tsv files
#' @param metadata Sample metadata. Path to location or object.
#' @param expr_mats Vector of expression matrices. Path to locations or objects.
#' @param covar_factors Vector of strings indicating additive covariates.
#' @param annots Annotations file. Path to location or object.
#' @param tissues Vector of strings indicating tissues in project. Ex: c("Kd","Lv") for "kidney" and "liver".
#' @param gridFile Genome Grid. Path to location or object. Defaults to 75k grid loaded with package.
#' @param suggLOD Suggestive LOD to use as filter for mediation. Default is 7.
#' @param localRange What is defined as "local". Default is 10e6.
#' @param total_cores Number of available cores to use. Default is NULL.
#' @param save_t Should files be saved, returned, or both. Default is "sr" (save and return). To save only use "so", to return only use "ro".
#'
#' @return A list containing \itemize{
#' \item{peaks_list}{Unfiltered peaks for each tissue.}
#' \item{maps_list}{List of objects associated with mapping. See [mapQTL] help for details.}
#' \item{res_list}{List containing mediation results for each tissue.}
#' \item{effects_res}{List of objects associated with effects. See [qtl_effects] help for details.}}
#' @export
#'
#'
#' @importFrom tibble lst
#'
#'
runQTL <- function(geno_out = "gbrs_interpolated_genoprobs.rds", peaks_out = "mm39_peaks.rds", map_out = "mm39_mapping.rds",
                   med_out = "mm39_mediation_res.rds", effects_out = "mm39_effects.rds", outdir, gbrs_fileLoc,
                   metadata, expr_mats, covar_factors, annots, tissues = c(),
                   gridFile = gridfile, suggLOD = 7, localRange = 10e6, total_cores = NULL, save_t = "sr") {
  ## Check oudir
  if (length(outdir) == 0 | !dir.exists(outdir)) {
    message("Invalid or no directory provided. Making an output file directory in the current working directory.")
    dir.create("./QTL_mapping")
    outdir <- "./QTL_mapping"
  }

  ## Check if object or file location is passed. If object assign to genoprobs, if file location run geoprobably
  if (is.list(gbrs_fileLoc)) {
    if (is.array(gbrs_fileLoc[[1]])) {
      message("converting to qtl2")
      genoprobs <- list()
      for (tissue in names(gbrs_fileLoc)) {
        genoprobs[[tissue]] <- probs_3d_to_qtl2(gbrs_fileLoc[[tissue]])
      }
    }
    if (is.list(gbrs_fileLoc[[1]])) {
      message("using provided probabilities")
      genoprobs <- gbrs_fileLoc
    }
  }
  if (is.character(gbrs_fileLoc)) {
    ## Convert genoprobs
    message("running genoprobs")
    genoprobs <- genoprobably(outfile = paste0(outdir, "/", geno_out), gbrsFileLoc = gbrs_fileLoc, tissues = tissues, gridFile = gridFile, save = save_t)
  }

  ## Map QTLs from genoprobs
  message("running mapping")
  map_peaks <- mapQTL(
    outdir = outdir, peaks_out = peaks_out, map_out = map_out, genoprobs = genoprobs,
    samp_meta = metadata, expr_mats = expr_mats, covar_factors = covar_factors,
    gridFile = gridFile, localRange = localRange, annots = annots, total_cores = total_cores, save = save_t
  )

  peaks_list <- map_peaks$peaks_list
  maps_list <- map_peaks$maps_list

  rm(map_peaks)

  ## Run Mediation and Effects
  message("running mediation")
  res_list <- run_mediate(peaks = peaks_list, mapping = maps_list, suggLOD = suggLOD, outdir = outdir, annots = annots, med_out = med_out, total_cores = total_cores, save = save_t)

  message("running effects")
  effects_res <- qtl_effects(mapping = maps_list, peaks = peaks_list, suggLOD = suggLOD, outdir = outdir, outfile = effects_out, total_cores = total_cores, save = save_t)

  if (save_t %in% c("sr", "ro")) {
    ## Return created objects
    all_out <- tibble::lst(peaks_list, maps_list, res_list, effects_res)
    return(all_out)
  }
}
