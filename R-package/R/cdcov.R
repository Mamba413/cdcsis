#' Conditional Distance Covariance with Given Bandwidth
#' 
#' Performs conditional distance covariance with given bandwidth.
#' 
#' Perform the measurement of the conditional independence between x and y
#' given z with given bandwidth.
#' 
#' @param x a numeric vector or matrix
#' @param y a numeric vector or matrix with compatible dimensions to \code{x}
#' @param z the variable being conditioned. \code{z} is a numeric vector or
#' matrix with compatible dimensions to \code{x}
#' @param width a positive value, user-specified
#' @param index exponent on Euclidean distance, in (0,2]
#' @return \item{cdcov}{conditional distance covariance with given variable
#' \code{z}; has the same length with \code{z}} \item{width}{the bandwidth}
#' \item{index}{ exponent on Euclidean distance, in (0,2] }
#' @author Canhong Wen, Wenliang Pan, and Xueqin Wang
#' @seealso \code{\link{cdcor}}
#' @references Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014).
#' Conditional distance correlation sure independence screening for ultrahigh
#' dimensional data. \emph{Submitted to Biostatistics}.
#' @keywords conditional distance covariance
#' @export
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- rnorm(100)
#' cdcov(x,y,z,0.25)
cdcov <- function(x, y, z, width, index = 1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  dim_x <- dim(x)
  n <- dim_x[1]
  p <- dim_x[2]
  q <- dim(y)[2]
  d <- dim(z)[2]
  k <- numeric(n * n)
  CDCOV <- numeric(n)
  re <-
    .C(
      "cdCOV",
      as.double(t(x)),
      as.double(t(y)),
      as.double(t(z)),
      as.integer(n),
      as.integer(p),
      as.integer(q),
      as.integer(d),
      as.double(index),
      as.double(width),
      k = as.double(k),
      cd = as.double(CDCOV)
    )
  cdc <- list(cdcov = re$cd,
              width = width,
              index = index)
  return(cdc)
}
