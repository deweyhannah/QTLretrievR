% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runQTL.R
\name{runQTL}
\alias{runQTL}
\title{Wrapper function to generate mapping, peaks, mediation, and effects data}
\usage{
runQTL(
  geno_out = "gbrs_interpolated_genoprobs.rds",
  peaks_out = "mm39_peaks.rds",
  map_out = "mm39_mapping.rds",
  med_out = "mm39_mediation_res.rds",
  effects_out = "mm39_effects.rds",
  outdir,
  gbrs_fileLoc,
  metadata,
  expr_mats,
  covar_factors,
  annots,
  tissues = c(),
  gridFile = gridfile,
  suggLOD = 7,
  localRange = 1e+07,
  total_cores = NULL,
  save_t = "sr"
)
}
\arguments{
\item{geno_out}{Output file name to save gbrs interpolated genoprobs. Default is "gbrs_interpolated_genoprobs.rds".}

\item{peaks_out}{Output file name to save peaks. Default is "mm39_peaks.rds".}

\item{map_out}{Output file name to save mapping. Default is "mm39_mapping.rds".}

\item{med_out}{Output file name to save mediation. Default is "mm39_mediation_res.rds".}

\item{effects_out}{Output file name to save effects. Default is "mm39_effects.rds".}

\item{outdir}{Output directory to save results.}

\item{gbrs_fileLoc}{Path to GBRS interpolated tsv files}

\item{metadata}{Sample metadata. Path to location or object.}

\item{expr_mats}{Vector of expression matrices. Path to locations or objects.}

\item{covar_factors}{Vector of strings indicating additive covariates.}

\item{annots}{Annotations file. Path to location or object.}

\item{tissues}{Vector of strings indicating tissues in project. Ex: c("Kd","Lv") for "kidney" and "liver".}

\item{gridFile}{Genome Grid. Path to location or object. Defaults to 75k grid loaded with package.}

\item{suggLOD}{Suggestive LOD to use as filter for mediation. Default is 7.}

\item{localRange}{What is defined as "local". Default is 10e6.}

\item{total_cores}{Number of available cores to use. Default is NULL.}

\item{save_t}{Should files be saved, returned, or both. Default is "sr" (save and return). To save only use "so", to return only use "ro".}
}
\value{
A list containing \itemize{
\item{peaks_list}{Unfiltered peaks for each tissue.}
\item{maps_list}{List of objects associated with mapping. See \link{mapQTL} help for details.}
\item{res_list}{List containing mediation results for each tissue.}
\item{effects_res}{List of objects associated with effects. See \link{qtl_effects} help for details.}}
}
\description{
Wrapper function to generate mapping, peaks, mediation, and effects data
}
