% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapQTL.R
\name{mapQTL}
\alias{mapQTL}
\title{Generate mapping data and peaks for QTL analysis}
\usage{
mapQTL(
  outdir,
  peaks_out,
  map_out,
  genoprobs,
  samp_meta,
  expr_mats,
  covar_factors,
  thrA = 5,
  thrX = 5,
  gridFile = gridfile,
  localRange = 1e+07,
  annots = NULL,
  total_cores = NULL,
  save = "sr"
)
}
\arguments{
\item{outdir}{Output directory where files mapping and peaks lists should be saved.}

\item{peaks_out}{String indicating the name for output peaks file. Should end in .rds.}

\item{map_out}{String indicating the name for output mapping file. Should end in .rds.}

\item{genoprobs}{Either a string with the name of the genoprobs file, or the genoprobs object.}

\item{samp_meta}{Sample metadata. Either a string pointing to the file, or the object itself.}

\item{expr_mats}{Vector of expression matrix files. One for each tissue, in the order that tissues were supplied to genoprobs.}

\item{covar_factors}{Additive covariate factors. These need to be columns in the factor metadata.}

\item{thrA}{Minimum reported LOD threshold for autosomes. Default is 5.}

\item{thrX}{Minimum reported LOD threshold for X chromosome. Default is 5.}

\item{gridFile}{File location for genome grid. Defaults to object loaded with package for 75k grid.}

\item{localRange}{What is defined as "local". Default is 10e6.}

\item{annots}{String pointing to annotations file or annotations object.}

\item{total_cores}{Number of available cores to use for parallelization. Default is NULL.}

\item{save}{Should files be saved, returned, or both. Default is "sr" (save and return). To save only use "so", to return only use "ro".}
}
\value{
A list containing: \itemize{
\item{maps_list}{A list of dataframes and lists to that can be used for future analyses and in other functions \itemize{
\item{qtlprobs}{Genome probabilities in qtl2format}
\item{covar_list}{list of covariate matrices for each tissue}
\item{expr_list}{Original normalized expression data for each tissue}
\item{exprZ_list}{Rank Z normalized expression data for each tissue}
\item{kinship_loco}{Kisnhip Matrix calculated using the "loco" option in \code{qtl2::calc_kinship}}
\item{gmap}{Genomic map of markers}
\item{map_dat2}{Combined genomic and physical map of markers}
\item{pmap}{Physical map of markers}
\item{tissue_samp}{Metadata broken down for each tissue}}}
\item{peaks_list}{A list of peaks list for each tissue.}}
}
\description{
Generate mapping data and peaks for QTL analysis
}
