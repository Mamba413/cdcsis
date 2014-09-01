\name{cdcor.ada}
\alias{cdcor.ada}
\title{
 Conditional Distance Correlation with Adaptive Bandwidth
}
\description{
Performs conditional distance correlation with adaptive bandwidth.
}
\usage{
cdcor.ada(x, y, z, tol = 0.1,index = 1)
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
  \item{tol}{
 the tolerance used in the bandwidth selection
}
  \item{index}{
 exponent on Euclidean distance, in (0,2] 
}
}
\details{
 It performs conditional distance correlation with adaptive bandwidth. The bandwidth is 
determined by the \code{bw} function.
}
\value{
  Returns an object with
  \item{cdcor }{conditional distance correlation}
  \item{width }{bandwidth used in \code{cor}. It is determined by the \code{bw} function.}

}
\references{
Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014). Conditional distance correlation 
sure independence screening for ultrahigh dimensional data. \emph{Submitted to Biostatistics}.
}
\author{
Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang
}
\seealso{
 \code{\link{cdcor}}
}
\examples{
 set.seed(1)
 # load the distance correlation for comparison
 require(energy) 

 ## independent case
 x <- rnorm(100)
 y <- rnorm(100)
 z <- rnorm(100)
 cdcor.ada(x,y,z)
 dcor(x,y)

 ## conditional dependent case
 x <- rnorm(100)
 y <- x + 0.1*rnorm(100)
 z <- rnorm(100)
 cdcor.ada(x,y,z)
 dcor(x,y)

 ## conditional independent case: x and z are correlated
 require(MASS)
 data <- mvrnorm(100,rep(0,2),matrix(c(1,0.8,0.8,1),2,2))
 x <- data[,1]
 z <- data[,2]
 y <- z + 0.1*rnorm(100)
 cdcor.ada(x,y,z)
 dcor(x,y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ conditional distance correlation }% __ONLY ONE__ keyword per line
