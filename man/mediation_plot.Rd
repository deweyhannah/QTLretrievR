% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediation_plot.R
\name{mediation_plot}
\alias{mediation_plot}
\title{Plot the change in mediation scores for the top 5}
\usage{
mediation_plot(
  effects,
  mediation,
  trans.bands,
  annots = annot_105,
  tissue,
  peak.chr,
  peak.num,
  outdir = NULL,
  sigLOD = 7.5,
  local = 5e+06,
  psave = FALSE
)
}
\arguments{
\item{effects}{Effects object from \code{qtl_effects()}}

\item{mediation}{Mediation object from \code{run_mediate()}}

\item{trans.bands}{Transbands object from \code{transbands()}}

\item{annots}{Molecular annotations. Default is ensembl v105 annotations.}

\item{tissue}{The tissue your transband of interest is in}

\item{peak.chr}{The chromosome your transband of interest is on}

\item{peak.num}{The transband number. If there is more than one transband on your chromosome of interest, specify which one (left to right) if there is only one then set this to 1}

\item{outdir}{Output directory where the plot will be saved. Only required if \code{psave = TRUE}}

\item{sigLOD}{Significant LOD threshold, determined by \code{LOD_thld()}}

\item{local}{How are you defining local in the transband. Default is 5e6}

\item{psave}{Should the plot be saved, or just returned? Default is FALSE}
}
\value{
A list containing: \itemize{
\item{results}{The mediation results from your hotspot with ranks}
\item{meds_ranked_sum}{A summary of the mediators and the number of times they are the top mediator for any given target}
\item{p}{Mediation Plot showing the LOD drop of the top 5 mediators for each target eQTL}}
}
\description{
Plot the change in mediation scores for the top 5
}
