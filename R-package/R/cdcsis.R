#' @title Conditional Distance Correlation Sure Independence Screening (CDC-SIS)
#' 
#' @description Performs conditional distance correlation sure independence screening (CDC-SIS).
#' 
#' @inheritParams cdcov.test
#' @param x a numeric matrix, or a list which contains multiple numeric matrix
#' @param distance if \code{distance = TRUE}, only \code{y} will be considered as distance matrices. Default: \code{distance = FALSE}
#' @param threshold the threshold of the number of predictors recuited by CDC-SIS. 
#' Should be less than or equal than the number of column of \code{x}. Default value \code{threshold} is sample size.
#' @return 
#' \item{ ix }{the vector of indices selected by CDC-SIS } 
#' \item{ cdcor }{ the conditional distance correlation for each univariate/multivariate variable in \code{x}} 
#' @author Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang
#' @seealso \code{\link{cdcor}}
#' @references Wen, C., Pan, W., Huang, M. and Wang, X., 2018. Sure independence screening adjusted for confounding covariates with ultrahigh-dimensional data. Statistica Sinica, 28, pp.293-317. URL http://www3.stat.sinica.edu.tw/statistica/J28N1/J28N114/J28N114.html
#' @export
#' @examples
#' \dontrun{
#' 
#' library(cdcsis)
#' 
#' ########## univariate explanative variables ##########
#' set.seed(1)
#' num <- 100
#' p <- 150
#' x <- matrix(rnorm(num * p), nrow = num)
#' z <- rnorm(num)
#' y <- 3 * x[, 1] + 1.5 * x[, 2] + 4 * z * x[, 5] + rnorm(num)
#' res <- cdcsis(x, y, z)
#' head(res[["ix"]], n = 10)
#' 
#' ########## multivariate explanative variables ##########
#' x <- as.list(as.data.frame(x))
#' x <- lapply(x, as.matrix)
#' x[[1]] <- cbind(x[[1]], x[[2]])
#' x[[2]] <- NULL
#' res <- cdcsis(x, y, z)
#' head(res[["ix"]], n = 10)
#' 
#' ########## multivariate response variables ##########
#' num <- 100
#' p <- 150
#' x <- matrix(rnorm(num * p), nrow = num)
#' z <- rnorm(num)
#' y1 <- 3 * x[, 1] + 5 * z * x[, 4] + rnorm(num)
#' y2 <- 3 * x[, 2] + 5 * x[, 3] + 2 * z + rnorm(num)
#' y <- cbind(y1, y2)
#' res <- cdcsis(x, y, z)
#' head(res[["ix"]], n = 10)
#' }
cdcsis <- function(x, y, z = NULL, 
                   kernel.type = c("gauss", "rectangle"), width, k,
                   threshold = nrow(y), hergenerity.z = FALSE, 
                   distance = FALSE, index = 1, num.threads = 1) 
{
  conditional.distance <- FALSE
  if (length(kernel.type) > 1) {
    kernel.type <- "gauss"
  } else {
    kernel.type <- match.arg(kernel.type)
  }
  z <- as.matrix(z)
  check_xyz_arguments(z)
  
  check_k_arguments(k, conditional.distance, kernel.type, z)
  
  # check heterogeneity and discrete
  if (!hergenerity.z) {
    hergenerity.z <- check_z_hergenerity(z)
  }
  discrete <- any_discrete_variable(z)
  
  if (any(discrete)) {
    if (hergenerity.z) {
      width <- discrete_z_heterogeneity_width(z, discrete)
    } else {
      width <- discrete_z_width(z)
    }
  } else {
    if (conditional.distance) {
      z <- compute_distance_matrix(z, conditional.distance, index)
      check_xyz_arguments(z)
      kernel.type <- "rectangle"
      width <- distance_z_width(z, k)
    } else {
      z <- as.matrix(z)
      check_xyz_arguments(z)
      if (kernel.type == "gauss") {
        if (missing(width)) {
          width <- gauss_z_width(z)
        }
      } else {
        if (missing(k)) {
          k <- compute_k(k, z)
        }
        if (hergenerity.z) {
          # width <- rectangle_z_hergenerity_width(z, k)
          width <- rectangle_z_hergenerity_width2(z)
        } else {
          # width <- rectangle_z_width(z, k) 
          width <- rectangle_z_width2(z) 
        }
      }
    } 
  }
  check_width_arguments(width)
  width <- as.double(width)
  
  check_threads_arguments(num.threads)
  check_index_arguments(index)
  
  y <- as.matrix(y)
  y <- compute_distance_matrix(y, distance, index)
  check_xyz_arguments(y)
  
  if (is.list(x)) {
    variable_index <- cumsum(sapply(x, ncol)) - 1
    names(variable_index) <- NULL
    variable_index <- as.integer(variable_index)
    x <- do.call("cbind", x)
    x <- as.matrix(x)
  } else {
    variable_index <- integer(0)
  }
  check_xyz_arguments(x)
  
  check_sample_size(y, z)
  check_sample_size(x, z)
  check_threshold_arguments(threshold, ncol(x))
  
  x <- as.matrix(t(x))
  
  kernel.type <- ifelse(kernel.type == "gauss", 1, 3)
  res <- cdcsisCpp(2, x, variable_index, y, z, width, index, num.threads, 0, 0, 2,
                   as.integer(kernel.type), as.integer(conditional.distance))
  res <- res[["statistic"]]
  list("ix" = order(res, decreasing = T)[1:threshold], 
       "cdcor" = res)
}