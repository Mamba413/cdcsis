#' @importFrom ks hpi
#' @importFrom ks Hpi.diag
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dist
#' @useDynLib cdcsis, .registration = TRUE
NULL


#' Convert to distance matrix
#'
#' @param x 
#' @param distance 
#' @param index 
#' @return distance matrix
#' @noRd
compute_distance_matrix <- function(x, distance, index) {
  x <- as.matrix(x)
  if (!distance) {
    x <- as.matrix(dist(x))
    if (index != 1) {
      x <- x^index
    }
  }
  x
}