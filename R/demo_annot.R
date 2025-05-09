#' Demo Annotations
#'
#' Ensembl v84 annotations to be used for vignette.
#'
#' @docType data
#' @name demo_annot
#' @format A 46521 x 6 data frame with gene annotations from Ensembl v84.
#' \describe{
#'   \item{id}{Ensembl gene IDs}
#'   \item{symbol}{Associated MGI symbol}
#'   \item{chr}{Chromosome of the gene}
#'   \item{start}{Start position (bp) on chromosome}
#'   \item{end}{End position (bp) on chromosome}
#'   \item{strand}{What strand the gene is on}
#' }
#'
#' @source Skelly et al., Mapping the Effects of Genetic Variation on Chromatin State and Gene Expression Reveals Loci That Control Ground State Pluripotency, Cell Stem Cell (2020), https://doi.org/10.1016/j.stem.2020.07.005
#' @usage data(demo_annot)
#'
NULL
