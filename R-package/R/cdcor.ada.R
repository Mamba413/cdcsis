#' Conditional Distance Correlation with Adaptive Bandwidth
#' 
#' Performs conditional distance correlation with adaptive bandwidth.
#' 
#' It performs conditional distance correlation with adaptive bandwidth. The
#' bandwidth is determined by the \code{bw} function.
#' 
#' @param x a numeric vector or matrix
#' @param y a numeric vector or matrix with compatible dimensions to \code{x}
#' @param z the variable being conditioned. \code{z} is a numeric vector or
#' matrix with compatible dimensions to \code{x}
#' @param tol the tolerance used in the bandwidth selection
#' @param index exponent on Euclidean distance, in (0,2]
#' @return Returns an object with \item{cdcor }{conditional distance
#' correlation} \item{width }{bandwidth used in \code{cor}. It is determined by
#' the \code{bw} function.}
#' @author Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang
#' @seealso \code{\link{cdcor}}
#' @references Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014).
#' Conditional distance correlation sure independence screening for ultrahigh
#' dimensional data. \emph{Submitted to Biostatistics}.
#' @keywords conditional distance correlation
#' @export
#' @examples
#' set.seed(1)
#' # load the distance correlation for comparison
#' require(energy) 
#' 
#' ## independent case
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- rnorm(100)
#' cdcor.ada(x,y,z)
#' dcor(x,y)
#' 
#' ## conditional dependent case
#' x <- rnorm(100)
#' y <- x + 0.1*rnorm(100)
#' z <- rnorm(100)
#' cdcor.ada(x,y,z)
#' dcor(x,y)
#' 
#' ## conditional independent case: x and z are correlated
#' require(MASS)
#' data <- mvrnorm(100,rep(0,2),matrix(c(1,0.8,0.8,1),2,2))
#' x <- data[,1]
#' z <- data[,2]
#' y <- z + 0.1*rnorm(100)
#' cdcor.ada(x,y,z)
#' dcor(x,y)
cdcor.ada <- function(x, y, z, tol = 1e-1, index = 1) {
  if (is.null(dim(z)) == TRUE) {
    width <- bw(x, y , z, index)
    if (width < tol)
      width <- hpi(z)
    out <- cdcor(x, y, z, width, index)$mcdcor
  }
  
  if (is.null(dim(z)) == FALSE) {
    width <- rep(NA, dim(z)[2])
    for (i in 1:dim(z)[2]) {
      width[i] <- bw(x, y, z[, i], index)
      if (width[i] < tol)
        width[i] <- hpi(z[, i])
    }
    w <- diag(width)
    out <- cdcor(x, y, z, w, index)$mcdcor
  }
  return(list(cdcor = out, width = width))
}
