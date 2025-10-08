#' Wrapper function to generate mapping, peaks, mediation, and effects data
#'
#' @param geno_out String indicating the name of the output file containing
#' interpolated genotype probabilities. This file will be saved in `.rds`
#' format and used for QTL mapping. Should end in `.rds`. Default is
#' "`gbrs_interpolated_genoprobs.rds`".
#' @param peaks_out String indicating the name for output peaks file.
#' This file will be saved in `.rds` format and will be used as an input
#' for downstream analysis. Should end  in `.rds`. Default is "`peaks.rds`"
#' @param map_out String indicating the name of the output file containing QTL
#' mapping results. This file will be saved in `.rds` format and will be used
#' in downstream analyses. Should end in `.rds`. Default is "`map.rds`"
#' @param med_out String indicating the name of the output file containing
#' mediation results for mediation within a phenotype. This file will be saved
#'  in `.rds` format and used for downstream analysis and visualization. Should
#'   end in `.rds`. Default is "`mediation.rds`"
#' @param effects_out String indicating the name of the output file containing
#'  founder haplotype effects results. This file will be saved in `.rds`
#'  format and used for downstream analysis and visualization. Should end in
#'  `.rds`. Default is "`effects.rds`"
#' @param outdir Directory to save output files. Default is `NULL`.
#' @param gbrs_fileLoc Path to GBRS interpolated tsv files.
#' @param metadata Sample metadata. Either a string pointing to the file,
#'  or the object itself.
#' @param expr_mats List of normalized count matrices (objects), or character
#'  paths to the file. One matrix per tissue. The order *must match* the
#'  tissue order in `genoprobs`.
#' @param covar_factors Additive covariate factors. These need to be columns
#'  in the sample metadata.
#' @param annots Annotations file. Contains mapping information for phenotypes.
#'  Dataframe, or tsv. Columns must include "id", "symbol", "start", "end".
#' @param tissues Vector of strings indicating tissues or conditions in project.
#' @param gridFile Genome Grid. Path to location or object. Defaults to
#' 75k grid loaded with package.
#' @param total_cores Number of available cores to use for parallelization.
#' Default is `NULL`.
#' @param save_t Indicates object return/save behavior. One of
#' `c("sr", "so", "ro")`; save & return, save only, return only.
#'  Default is "sr".
#'
#' @return A list containing \itemize{
#' \item{peaks_list}{Unfiltered peaks for each tissue.}
#' \item{maps_list}{List of objects associated with mapping. See [mapQTL]
#' help for details.}
#' \item{res_list}{List containing mediation results for each tissue.}
#' \item{effects_res}{List of objects associated with effects. See [qtl_effects]
#' help for details.}}
#' @export
#'
#'
#' @importFrom tibble lst
#'
#'
runQTL <- function(geno_out = "gbrs_interpolated_genoprobs.rds",
                   peaks_out = "peaks.rds", map_out = "mapping.rds",
                   med_out = "mediation_res.rds",
                   effects_out = "effects.rds", outdir = NULL, gbrs_fileLoc,
                   metadata, expr_mats, covar_factors, annots, tissues = c(),
                   gridFile = gridfile, total_cores = NULL,
                   save_t = "sr", ...) {
  ## Check outdir
  if(save_t %in% c("sr", "so")) {
    if (is.null(outdir)) {
      stop("Files to be saved, but no output directory provided. Please
           either provide a valid path for outdir, or change save_t to 'ro'")
    }
    if (!is.null(outdir) & !dir.exists(outdir)) {
      message("Invalid directory provided. Making an output file directory
            in the current working directory.")
      dir.create("./QTL_mapping")
      outdir <- "./QTL_mapping"
    }
  }

  save_og <- save_t

  if (save_og == "so") {
    save_t <- "sr"
  }

  opt_args <- list(...)
  if("thrA" %notin% names(opt_args)) {
    thrA <- 5
  } else {
    thrA <- opt_args$thrA
  }
  if("thrX" %notin% names(opt_args)) {
    thrX <- 5
  } else {
    thrX <- opt_args$thrX
  }
  if("localRange" %notin% names(opt_args)) {
    localRange <- 2e6
  } else {
    localRange <- opt_args$localRange
  }
  if("rz" %notin% names(opt_args)) {
    rz <- FALSE
  } else {
    rz <- opt_args$rz
  }
  if("phys" %notin% names(opt_args)) {
    phys <- TRUE
  } else {
    phys <- opt_args$phys
  }
  ## Check if object or file location is passed. If object assign to genoprobs,
  ## if file location run geoprobably
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
    genoprobs <- genoprobably(outfile = paste0(outdir, "/", geno_out),
                              gbrsFileLoc = gbrs_fileLoc, tissues = tissues,
                              gridFile = gridFile, save = save_t)
  }

  ## Map QTLs from genoprobs
  message("running mapping")
  map_peaks <- mapQTL(
    outdir        = outdir,
    peaks_out     = peaks_out,
    map_out       = map_out,
    genoprobs     = genoprobs,
    samp_meta     = metadata,
    expr_mats     = expr_mats,
    covar_factors = covar_factors,
    gridFile      = gridFile,
    annots        = annots,
    total_cores   = total_cores,
    save          = save_t,
    thrA          = thrA,
    thrX          = thrX,
    localRange    = localRange,
    rz            = rz,
    phys          = phys,
    ...
    )

  peaks_list <- map_peaks$peaks_list
  maps_list <- map_peaks$maps_list

  rm(map_peaks)

  ## Run Mediation and Effects
  message("running mediation")
  res_list <- modiFinder(peaks       = peaks_list,
                         mapping     = maps_list,
                         sigLOD      = sigLOD,
                         outdir      = outdir,
                         annots      = annots,
                         med_out     = med_out,
                         total_cores = total_cores,
                         save        = save_t)

  message("running effects")
  effects_res <- qtl_effects(mapping         = maps_list,
                             peaks           = peaks_list,
                             suggLOD         = sigLOD,
                             outdir          = outdir,
                             effects_out     = effects_out,
                             total_cores     = total_cores,
                             save            = save_t)

  save_t <- save_og

  if (save_t %in% c("sr", "ro")) {
    ## Return created objects
    all_out <- tibble::lst(peaks_list, maps_list, res_list, effects_res)
    return(all_out)
  }
}
