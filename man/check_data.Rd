% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_data.R
\name{check_data}
\alias{check_data}
\title{Checking the inputted data for QTLretrievR}
\usage{
check_data(x, type = "")
}
\arguments{
\item{x}{can be a path to an Rds object or a list of R objects.}

\item{type}{The name of the object if a single one is passed. Use "genoprobs" for genotype probabilities, "mediation" for mediation results and "peaks" for a table of eQTL peaks.}
}
\value{
List of properly names R objects to be used in downstream analysis of QTLretrievR.
}
\description{
Checking the inputted data for QTLretrievR
}
