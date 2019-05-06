#' @title Conditional Distance Covariance/Correlation Statistics
#' @description Computes conditional distance covariance and conditional distance correlation statistics, 
#' which are multivariate measures of conditional dependence.
#' 
#' @inheritParams cdcov.test
#' @param width a user-specified positive value (univariate conditional variable) or vector (multivariate conditional variable) for 
#' gaussian kernel bandwidth. Its default value is relies on \code{stats::bw.nrd0}.
#' @rdname cdcov
#' @details 
#' \code{cdcov} and \code{cdcor} compute conditional distance covariance and conditional distance correlation statistics.
#' The sample sizes (number of rows or length of the vector) of the two variables must agree, 
#' and samples must not contain missing values. 
#' If we set \code{distance = TRUE}, arguments \code{x}, \code{y} can be a \code{dist} object recording distance between samples; 
#' otherwise, these arguments are treated as multivariate data.
#' 
#' @return 
#' \item{\code{cdcov }}{ sample version of conditional distance covariance.}
#' 
#' @author Canhong Wen, Wenliang Pan, and Xueqin Wang
#' @seealso \code{\link{cdcor}}
#' @references Wang, X., Pan, W., Hu, W., Tian, Y. and Zhang, H., 2015. Conditional distance correlation. Journal of the American Statistical Association, 110(512), pp.1726-1734.
#' 
#' @keywords conditional distance covariance
#' @export
#' 
#' @examples
#' library(cdcsis)
#' 
#' ############# Conditional Distance Correlation #############
#' set.seed(1)
#' x <- rnorm(25)
#' y <- rnorm(25)
#' z <- rnorm(25)
#' cdcov(x, y, z)
#' 
cdcov <- function(x, y, z, kernel.type = c("rectangle", "gauss"), k = 6, width,
                  index = 1, distance = FALSE) 
{
  conditional.distance <- FALSE
  if (length(kernel.type) > 1) {
    kernel.type <- "gauss"
  } else {
    kernel.type <- match.arg(kernel.type)
  }
  
  z <- as.matrix(z)
  check_xyz_arguments(z)
  if (missing(width)) {
    if (dim(z)[2] == 1) {
      width <- stats::bw.nrd0(as.vector(z))
    } else if (dim(z)[2] <= 3) {
      width <- diag(ks::Hpi.diag(z))
    } else {
      width <- apply(z, 2, stats::bw.nrd)
    }
  }
  check_width_arguments(width)
  width <- as.double(width)
  
  check_index_arguments(index)
  
  x <- compute_distance_matrix(x, distance, index)
  check_xyz_arguments(x)
  
  y <- compute_distance_matrix(y, distance, index)
  check_xyz_arguments(y)
  
  check_sample_size(x, z)
  check_sample_size(y, z)
  
  kernel.type <- ifelse(kernel.type == "gauss", 1, 3)
  res <- cdcsisCpp(stats_method = 3, x, c(0), y, z, width, index, 0, 0, 0, 1, 
                   as.integer(kernel.type), as.integer(conditional.distance))
  res <- res[["statistic"]]
  names(res) <- "cdcov"
  res
}

#' @inheritParams cdcov.test
#' @rdname cdcov
#' @return 
#' \item{\code{cdcor }}{ sample version of conditional distance correlation.}
#' @export
#' 
#' @examples
#' ############# Conditional Distance Correlation #############
#' num <- 25
#' set.seed(1)
#' x <- rnorm(num)
#' y <- rnorm(num)
#' z <- rnorm(num)
#' cdcor(x, y, z)
cdcor <- function(x, y, z, kernel.type = c("rectangle", "gauss"), k = 6, width, 
                  index = 1, distance = FALSE) {
  conditional.distance <- FALSE
  if (length(kernel.type) > 1) {
    kernel.type <- "gauss"
  } else {
    kernel.type <- match.arg(kernel.type)
  }
  
  z <- as.matrix(z)
  check_xyz_arguments(z)
  if (missing(width)) {
    if (dim(z)[2] == 1) {
      width <- stats::bw.nrd0(as.vector(z))
    } else if (dim(z)[2] <= 3) {
      width <- diag(ks::Hpi.diag(z))
    } else {
      width <- apply(z, 2, stats::bw.nrd)
    }
  }
  check_width_arguments(width)
  width <- as.double(width)
  
  check_index_arguments(index)
  
  x <- compute_distance_matrix(x, distance, index)
  check_xyz_arguments(x)
  
  y <- compute_distance_matrix(y, distance, index)
  check_xyz_arguments(y)
  
  check_sample_size(x, z)
  check_sample_size(y, z)
  
  kernel.type <- ifelse(kernel.type == "gauss", 1, 3)
  res <- cdcsisCpp(stats_method = 3, x, c(0), y, z, width, index, 0, 0, 0, 2,
                   as.integer(kernel.type), as.integer(conditional.distance))
  res <- res[["statistic"]]
  names(res) <- "cdcor"
  res
}


#' @title Conditional Distance Covariance Independence Test
#' 
#' @description Performs the nonparametric conditional distance covariance test for conditional independence assumption
#'
#' @param x a numeric vector, matrix, or \code{dist} object
#' @param y a numeric vector, matrix, or \code{dist} object
#' @param z \code{z} is a numeric vector or matrix. It is the variable being conditioned.
#' @param num.bootstrap the number of local bootstrap procedure replications. Default: \code{num.bootstrap = 99}.
#' @param kernel.type the kernel to be used. This must be one of "gauss", "rectangle". Any unambiguous substring can be given.
#' Default: \code{kernel.type = "gauss"}.
#' @param width a user-specified positive value (univariate conditional variable) or vector (multivariate conditional variable) for 
#' gaussian kernel bandwidth. Its default value is relies on \code{stats::bw.nrd0} function when conditional variable is univariate, 
#' \code{ks::Hpi.diag} when conditional variable with at most trivariate, and \code{stats::bw.nrd} on the other cases.
#' @param k a user-specified positive integer value for computing rectangle kernel bandwidth.
#' @param index exponent on Euclidean distance, in \eqn{(0,2]}
#' @param distance if \code{distance = TRUE}, \code{x} and \code{y} will be considered as distance matrices. Default: \code{distance = FALSE}.
#' @param seed the random seed
#' @param num.threads number of threads. Default \code{num.threads = 1}.
#'
#' @return \code{cdcov.test} returns a list with class "htest" containing the following components:
#' \item{\code{statistic}}{conditional distance covariance statistic.}            
#' \item{\code{p.value}}{the \eqn{p}-value for the test.}
#' \item{\code{replicates}}{the number of local bootstrap procedure replications.}
#' \item{\code{size}}{sample sizes.}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' 
#' @rdname cdcov.test
#' @seealso \code{\link{cdcov}}
#' 
#' @export
#' 
#' @examples
#' library(cdcsis)
#' set.seed(1)
#' num <- 50
#' ################# Conditional Independent #################
#' ## Case 1:
#' cov_mat <- matrix(c(1, 0.36, 0.6, 0.36, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
#' dat <- mvtnorm::rmvnorm(n = num, sigma = cov_mat)
#' x <- dat[, 1]
#' y <- dat[, 2]
#' z <- dat[, 3]
#' cdcov.test(x, y, z)
#' ## Case 2:
#' z <- rnorm(num)
#' x <- 0.5 * (z^3 / 7 + z / 2) + tanh(rnorm(num))
#' x <- x + x^3 / 3
#' y <- (z^3 + z) / 3 + rnorm(num)
#' y <- y + tanh(y / 3)
#' cdcov.test(x, y, z, num.bootstrap = 99)
#' 
#' ################# Conditional Dependent #################
#' ## Case 3:
#' cov_mat <- matrix(c(1, 0.7, 0.6, 0.7, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
#' dat <- mvtnorm::rmvnorm(n = num, sigma = cov_mat)
#' x <- dat[, 1]
#' y <- dat[, 2]
#' z <- dat[, 3]
#' cdcov.test(x, y, z, width = 0.5)
#' ## Case 4:
#' z <- matrix(rt(num * 4, df = 2), nrow = num)
#' x <- z
#' y <- cbind(sin(z[, 1]) + cos(z[, 2]) + (z[, 3])^2 + (z[, 4])^2, 
#'            (z[, 1])^2 + (z[, 2])^2 + z[, 3] + z[, 4])
#' z <- z[, 1:2]
#' cdcov.test(x, y, z, seed = 2)
#' 
#' ################# Distance Matrix Input #################
#' x <- dist(x)
#' y <- dist(y)
#' cdcov.test(x, y, z, seed = 2, distance = TRUE)
cdcov.test <- function(x, y, z, num.bootstrap = 99, 
                       kernel.type = c("rectangle", "gauss"), k = 6, width,
                       distance = FALSE, index = 1, seed = 1, num.threads = 1) {
  conditional.distance <- FALSE
  if (length(kernel.type) > 1) {
    kernel.type <- "gauss"
  } else {
    kernel.type <- match.arg(kernel.type)
  }
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)), "and", deparse(substitute(z)))
  
  z <- as.matrix(z)
  check_xyz_arguments(z)
  
  if (missing(width)) {
    if (dim(z)[2] == 1) {
      width <- stats::bw.nrd0(as.vector(z))
    } else if (dim(z)[2] <= 3) {
      width <- diag(ks::Hpi.diag(z))
    } else {
      width <- apply(z, 2, stats::bw.nrd)
    }
  }
  check_width_arguments(width)
  width <- as.double(width)  
  x <- compute_distance_matrix(x, distance, index)
  check_xyz_arguments(x)
  
  y <- compute_distance_matrix(y, distance, index)
  check_xyz_arguments(y)
  
  check_sample_size(x, z)
  check_sample_size(y, z)
  
  kernel.type <- ifelse(kernel.type == "gauss", 1, 3)
  res <- cdcsisCpp(stats_method = 1, x, c(0), y, z, width, index, num.threads, num.bootstrap, seed, 1, 
                   kernel.type, as.integer(conditional.distance))
  
  res <- wrap_to_htest(res, num.bootstrap, nrow(z), data_name)
  names(res[["statistic"]]) <- "cdcov"
  res
}

#' @title Conditional Ball Covariance Test
#' @description Performs the nonparametric conditional ball covariance test for conditional independence assumption
#' @inheritParams cdcov.test
#' @param kernel.type the kernel to be used. This must be one of "gauss", "rectangle". Any unambiguous substring can be given.
#' Default: \code{kernel.type = "rectangle"}.
#' @rdname cbcov
#'
#' @return \code{cbcov.test} returns a list with class "htest" containing the following components:
#' \item{\code{statistic}}{conditional ball covariance statistic.}            
#' \item{\code{p.value}}{the \eqn{p}-value for the test.}
#' \item{\code{replicates}}{the number of local bootstrap procedure replications.}
#' \item{\code{size}}{sample sizes.}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' @export
#'
#' @examples
#' library(cdcsis)
#' set.seed(1)
#' num <- 50
#' ################# Conditional Independent #################
#' ## Case 1:
#' cov_mat <- matrix(c(1, 0.36, 0.6, 0.36, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
#' dat <- mvtnorm::rmvnorm(n = num, sigma = cov_mat)
#' x <- dat[, 1]
#' y <- dat[, 2]
#' z <- dat[, 3]
#' cbcov.test(x, y, z)
#' cbcov.test(x, y, z, kernel.type = "gauss")
#' ## Case 2:
#' z <- rnorm(num)
#' x <- 0.5 * (z^3 / 7 + z / 2) + tanh(rnorm(num))
#' x <- x + x^3 / 3
#' y <- (z^3 + z) / 3 + rnorm(num)
#' y <- y + tanh(y / 3)
#' cbcov.test(x, y, z, num.bootstrap = 199)
#' 
#' 
#' ################# Conditional Dependent #################
#' ## Case 3:
#' cov_mat <- matrix(c(1, 0.7, 0.6, 0.7, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
#' dat <- mvtnorm::rmvnorm(n = num, sigma = cov_mat)
#' x <- dat[, 1]
#' y <- dat[, 2]
#' z <- dat[, 3]
#' cbcov.test(x, y, z, k = 10)
#' ## Case 4:
#' z <- matrix(rt(num * 4, df = 2), nrow = num)
#' x <- z
#' y <- cbind(sin(z[, 1]) + cos(z[, 2]) + (z[, 3])^2 + (z[, 4])^2, 
#'            (z[, 1])^2 + (z[, 2])^2 + z[, 3] + z[, 4])
#' z <- z[, 1:2]
#' cbcov.test(x, y, z, seed = 2, k = 10)
#' 
#' ################# Distance Matrix Input #################
#' x <- dist(x)
#' y <- dist(y)
#' cbcov.test(x, y, z, seed = 1, distance = TRUE, k = 10)
cbcov.test <- function(x, y, z, num.bootstrap = 99, 
                       kernel.type = c("rectangle", "gauss"), k = 6, width,
                       index = 1, distance = FALSE, seed = 1, num.threads = 1) {
  conditional.distance <- FALSE
  if (length(kernel.type) > 1) {
    kernel.type <- "rectangle"
  } else {
    kernel.type <- match.arg(kernel.type)
  }
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)), "and", deparse(substitute(z)))
  
  if (conditional.distance) {
    z <- compute_distance_matrix(z, conditional.distance, index)
    check_xyz_arguments(z)
    width <- t(apply(z, 1, sort))
    width <- width[, k]
    kernel.type <- "rectangle"
  } else {
    z <- as.matrix(z)
    check_xyz_arguments(z)
    if (kernel.type == "gauss") {
      if (missing(width)) {
        if (dim(z)[2] == 1) {
          width <- stats::bw.nrd0(as.vector(z))
        } else if (dim(z)[2] <= 3) {
          width <- diag(ks::Hpi.diag(z))
        } else {
          width <- apply(z, 2, stats::bw.nrd)
        }
      }
    } else {
      width <- as.matrix(dist(z))
      width <- t(apply(width, 1, sort))
      width <- mean(width[, k]) / sqrt(dim(z)[2])
    }
  }
  
  check_width_arguments(width)
  width <- as.double(width)
  
  x <- compute_distance_matrix(x, distance, index)
  check_xyz_arguments(x)
  
  y <- compute_distance_matrix(y, distance, index)
  check_xyz_arguments(y)
  
  check_sample_size(x, z)
  check_sample_size(y, z)
  
  set.seed(seed)
  kernel.type <- ifelse(kernel.type == "gauss", 1, 3)
  res <- cdcsisCpp(stats_method = 1, x, c(0), y, z, width, index, num.threads, 
                   num.bootstrap, seed, 3, kernel.type, as.integer(conditional.distance))
  
  res <- wrap_to_htest(res, num.bootstrap, nrow(z), data_name, cdc = FALSE)
  names(res[["statistic"]]) <- "cbcov"
  res
}


wrap_to_htest <- function(res, num.bootstrap, num, data_name, cdc = TRUE) {
  res[["replicates"]] <- num.bootstrap
  res[["size"]] <- num
  res[["alternative"]] <- "random variables are conditional dependent"
  res[["method"]] <- sprintf("Conditional %s Covariance Test", ifelse(cdc, "Distance", "Ball"))
  data_name <- paste(data_name, sprintf("\nnumber of observations = %s, ", num))
  data_name <- paste0(data_name, "replicates = ", num.bootstrap)
  res[["data.name"]] <- data_name
  class(res) <- "htest"
  res
}
