#' Conditional Distance Correlation Sure Independence Screening (CDCSIS)
#' 
#' Performs conditional distance correlation sure independence screening
#' (CDCSIS).
#' 
#' It performs conditional distance correlation sure independence screening
#' (CDCSIS).
#' 
#' @param x a matrix
#' @param y a numeric vector or matrix with compatible dimensions to \code{x}
#' @param z the variable being conditioned. \code{z} is a numeric vector or
#' matrix with compatible dimensions to \code{x}
#' @param thres the threshold of the number of pedictors recuited by CDCSIS.
#' Should be less than or equal than the number of column of \code{x}.
#' @return \item{ CDCSISind }{the vector of indices selected by CDCSIS } \item{
#' thres }{the threshold of the number of pedictors recuited by CDCSIS} \item{
#' DC }{ the distance correlation for each dimensionality of \code{x}} \item{
#' DCord }{the order of \code{DC} for each dimensionality of \code{x} }
#' @author Canhong Wen, Wenliang Pan, and Xueqin Wang
#' @seealso \code{\link{cdcor.ada}}, \code{\link{cdcor}}
#' @references Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014).
#' Conditional distance correlation sure independence screening for ultrahigh
#' dimensional data. \emph{Submitted to Biostatistics}.
#' @keywords conditional distance correlation sure independence screening
#' Ultrahigh dimensional
#' @export
#' @examples
#' set.seed(0)
#' n <- 100
#' p <- 10
#' rho <- 0.5
#' Sigma <- matrix(rho, p, p)
#' diag(Sigma) <- 1
#' 
#' require(MASS)
#' x <- mvrnorm(n,rep(0,p), Sigma)
#' y <- x[,2] + x[,3] + rnorm(n)
#' z <- x[,1]
#' 
#' cdcsis(x,y,z,2)
#' 
cdcsis <- function(x, y, z, thres) {
  p <- dim(x)[2]
  DW.c <- sapply(1:p, function(j)
    cdcor.ada(x[, j], y, z))
  DCW <- unlist(DW.c[1, ])
  DCWord <- order(abs(DCW), decreasing = T)
  CDCSISind <- DCWord[1:thres]
  
  return(list(
    CDCSISind = CDCSISind,
    thres = thres,
    DC = DCW,
    DCord = DCWord
  ))
}


#' @import ks
#' @useDynLib cdcsis
NULL