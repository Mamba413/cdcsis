\name{cdcor}
\alias{cdcor}
\title{
 Conditional Distance Correlation with Given Bandwidth
}
\description{
Performs conditional distance correlation with given bandwidth.
}
\usage{
cdcor(x, y, z, width, index = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a numeric vector or matrix
}
  \item{y}{
 a numeric vector or matrix with compatible dimensions to \code{x}
}
  \item{z}{
 the variable being conditioned. \code{z} is a numeric vector or matrix 
 with compatible dimensions to \code{x}
}
  \item{width}{
 a positive value, user-specified
}
  \item{index}{
 exponent on Euclidean distance, in (0,2] 
}
}
\details{
It performs conditional distance correlation with given bandwidth.
}
\value{
   \item{ cdcov }{the conditional distance covariance with given variable z } 
   \item{ mcdcov }{mean of the conditional distance covariance with given variable z} 
    \item{ width }{ the bandwidth }
    \item{ index }{ exponent on Euclidean distance, in (0,2] }
}
\references{
Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014). Conditional distance correlation 
sure independence screening for ultrahigh dimensional data. \emph{Submitted to Biostatistics}.
}
\author{
Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang
}

\seealso{
 \code{\link{cdcor.ada}}, \code{\link{cdcov}}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
 x <- rnorm(100)
 y <- rnorm(100)
 z <- rnorm(100)
 cdcor(x,y,z,0.25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ conditional distance correlation }

