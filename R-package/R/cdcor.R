#' @title Conditional Distance Correlation with Given Bandwidth
#' @description Performs conditional distance correlation with given bandwidth. It performs conditional distance correlation with given bandwidth.
#' 
#' @param x a numeric vector or matrix
#' @param y a numeric vector or matrix with compatible dimensions to \code{x}
#' @param z the variable being conditioned. \code{z} is a numeric vector or matrix with compatible dimensions to \code{x}
#' @param width a positive value, user-specified
#' @param index exponent on Euclidean distance, in \eqn{(0,2]}
#' @return 
#' \item{cdcov}{the conditional distance covariance with given variable \code{z}}
#' \item{mcdcov}{mean of the conditional distance covariance with given variable \code{z}}
#' \item{width}{the bandwidth}
#' \item{index}{exponent on Euclidean distance, in \eqn{(0,2]}}
#' @author Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang
#' @seealso \code{\link{cdcor.ada}}, \code{\link{cdcov}}
#' @references Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014).
#' Conditional distance correlation sure independence screening for ultrahigh
#' dimensional data. \emph{Submitted to Biostatistics}.
#' @keywords conditional distance correlation
#' @export
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- rnorm(100)
#' cdcor(x, y, z, 0.25)
#' # 
#' z <- cbind(rnorm(100), rnorm(100))
#' cdcor(x, y, z, diag(1, 2))
cdcor <- function(x,y,z,width,index=1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  p <- 1
  dim_x <- dim(x)
  n <- dim_x[1]
  M <- dim_x[2]
  q <- dim(y)[2]
  d <- dim(z)[2]
  k <- numeric(n * n)
  CDCOV <- numeric(n)
  iraCDCOV <- numeric(M)
  re <-
    .C(
      "IterationcdCov",
      as.double((x)),
      as.double(t(y)),
      as.double(t(z)),
      as.integer(n),
      as.integer(p),
      as.integer(q),
      as.integer(d),
      as.double(index),
      as.double(width),
      k = as.double(k),
      cd = as.double(CDCOV),
      as.integer(M),
      iracd = as.double(iraCDCOV)
    )
  cdc <- list("cdcor" = re$cd, "mcdcor" = re$iracd,
      "width" = width, "index" = index)
  return(cdc)
}
