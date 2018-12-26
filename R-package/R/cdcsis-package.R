#' Conditional Distance Correlation and Its Related Feature Screening Method
#' 
#' Gives conditional distance correlation and performs the conditional distance correlation sure independence screening procedure for ultrahigh dimensional data. 
#' The conditional distance correlation is a novel conditional dependence measurement of two random variables given a third variable. 
#' The conditional distance correlation sure independence screening is used for screening variables in ultrahigh dimensional setting.
#' 
#' @name cdcsis-package
#' @docType package
#' @author Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang
#' 
#' Maintainer: Canhong Wen <wencanhong@@gmail.com>
#' @references Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014). Conditional distance correlation sure independence screening for ultrahigh dimensional data.
#' @keywords SIS, distance correlation, confounder
#' @examples
#' set.seed(0)
#' n <- 100 # sample size
#' p <- 10 # dimensionality
#' rho <- 0.5 # the correlation between pairwise predictors.
#' Sigma <- matrix(rho, p, p)
#' diag(Sigma) <- 1
#' 
#' require(MASS)
#' x <- mvrnorm(n,rep(0,p), Sigma)
#' 
#' y <- x[,2] + x[,3] + rnorm(n)
#' z <- x[,1]
#' cdcsis(x,y,z,2)
NULL