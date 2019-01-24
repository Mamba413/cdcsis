#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom KernSmooth dpik
#' @importFrom ks Hpi.diag
#' @importFrom stats dist
#' @importFrom stats bw.nrd0
#' @importFrom stats bw.nrd
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

#' check width arguments
#'
#' @param width 
#' @noRd
check_width_arguments <- function(width) {
  if (any(width <= 0)) {
    stop("width invalid is invalid!")
  }
}

check_xyz_arguments <- function(x) {
  if (anyNA(x)) {
    stop("input dataset include NA!")
  }
  if (any(is.infinite(x))) {
    stop("input dataset include infinite value!")
  }
}

#' check num.threads arguments
#' @param num.threads A integer number
#' @noRd
check_threads_arguments <- function(num.threads) {
  if(is.null(num.threads) | (num.threads < 1)) {
    stop("num.threads arguments is invalid!")
  }
}


#' check sample size of two dataset
#'
#' @param x numeric matrix
#' @param y numeric matrix
#' @noRd
check_sample_size <- function(x, y) {
  num_row_x <- dim(x)[1]
  num_row_y <- dim(y)[1]
  if(num_row_x != num_row_y) {
    stop("dataset with different sample size!")
  }
}


#' check num.threads arguments
#' @param num.threads A integer number
#' @noRd
check_index_arguments <- function(index) {
  if(is.null(index) | (index <= 0) | (index > 2)) {
    stop("index arguments is invalid!")
  }
}

#' check num.threads arguments
#' @param num.threads A integer number
#' @noRd
check_threshold_arguments <- function(threshold, num_col) {
  if(is.null(threshold) | (threshold <= 0) | (threshold > num_col)) {
    stop("threshold arguments is invalid!")
  }
}