#' Upper Quartile Normalization
#'
#' @param mat A matrix of raw counts to be normalized
#'
#' @return Upper Quartile Nomalized matrix of counts
#' @export
#'
uqn <- function(mat) {
  ## Calculate the 75th percentile
  if (!is.matrix(mat)) {
    stop("Input needs to be a matrix")
  }
  upper_quartiles <- apply(mat, 2, function(x) quantile(x[x > 0], 0.75, na.rm = TRUE))

  ## Divide by the upper quartiles
  normalized_counts <- sweep(mat, 2, upper_quartiles, FUN = "/")

  ## Fix the column and rownames
  rownames(normalized_counts) <- rownames(mat)
  colnames(normalized_counts) <- colnames(mat)

  return(normalized_counts)
}
