\name{cdcsis}
\alias{cdcsis}
\title{
 Conditional Distance Correlation Sure Independence Screening (CDCSIS)
}
\description{
Performs conditional distance correlation sure independence screening (CDCSIS).
}
\usage{
cdcsis(x, y, z, thres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a matrix
}
  \item{y}{
 a numeric vector or matrix with compatible dimensions to \code{x}
}
  \item{z}{
 the variable being conditioned. \code{z} is a numeric vector or matrix 
 with compatible dimensions to \code{x}
}
  \item{thres}{
the threshold of the number of pedictors recuited by CDCSIS. Should be less than 
or equal than the number of column of \code{x}.
}
}
\details{
It performs conditional distance correlation sure independence screening (CDCSIS).
}
\value{
   \item{ CDCSISind }{the vector of indices selected by CDCSIS } 
   \item{ thres }{the threshold of the number of pedictors recuited by CDCSIS} 
    \item{ DC }{ the distance correlation for each dimensionality of \code{x}}
    \item{ DCord }{the order of \code{DC} for each dimensionality of \code{x} }
}
\references{
Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014). Conditional distance correlation 
sure independence screening for ultrahigh dimensional data. \emph{Submitted to Biostatistics}.
}
\author{
Canhong Wen, Wenliang Pan, and Xueqin Wang
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{cdcor.ada}}, \code{\link{cdcor}}
}

\examples{
 set.seed(0)
 n <- 100
 p <- 10
 rho <- 0.5
 Sigma <- matrix(rho, p, p)
 diag(Sigma) <- 1

 require(MASS)
 x <- mvrnorm(n,rep(0,p), Sigma)
 y <- x[,2] + x[,3] + rnorm(n)
 z <- x[,1]

 cdcsis(x,y,z,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conditional distance correlation}
\keyword{ sure independence screening }
\keyword{ Ultrahigh dimensional }

