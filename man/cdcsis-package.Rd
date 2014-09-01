\name{cdcsis-package}
\alias{cdcsis-package}
\docType{package}
\title{
Conditional Distance Correlation and Its Related Feature Screening Method
}
\description{
Gives the conditional distance correlation and performs its based sure independence screening method, 
i.e., CDCSIS of Wen et al.(2014).
}
\details{
\tabular{ll}{
Package: \tab CDCSIS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-01\cr
License: \tab GPL(>=2) \cr
}
 The package \code{cdcsis} is used for calculate the conditional distance correlation and performs
 the related sure independent screening method. Details of the method can be found in Wen et.al. (2014).
}
\author{
Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang

Maintainer: Canhong Wen <wencanhong@gmail.com>
}
\references{
Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014). Conditional distance correlation 
sure independence screening for ultrahigh dimensional data. \emph{Submitted to Biostatistics}.
}
\keyword{ SIS, distance correlation, confounder }
\examples{
 set.seed(0)
 n <- 100 # sample size
 p <- 10 # dimensionality
 rho <- 0.5 # the correlation between pairwise predictors.
 Sigma <- matrix(rho, p, p)
 diag(Sigma) <- 1

 require(MASS)
 x <- mvrnorm(n,rep(0,p), Sigma)

 y <- x[,2] + x[,3] + rnorm(n)
 z <- x[,1]
 cdcsis(x,y,z,2)
}
