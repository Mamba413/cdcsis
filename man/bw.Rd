\name{bw}
\alias{bw}
\title{
  Bandwidth Selection in the Conditional Distance Correlation
}
\description{
\code{bw} is used to select the bandwidth in the conditional distance correlation estimation.
}
\usage{
bw(x, y, z, index = 1)
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
  \item{index}{
 exponent on Euclidean distance, in (0,2] 
}
}
\details{
 For univariate \code{z}, the univariate plug-in selector of Wand & Jones (1994) is used. That is,
 \code{hpi} function in the \code{ks} package is used.


 For multivariate \code{z}, an diagonal matrix of the bandwidth is assumed, i.e., select the 
 optimal bandwidth with \code{hpi} individually for each column of \code{z}.
}
\value{
  The plug-in bandwidth.
}
\references{
Canhong Wen, Wenliang Pan, Mian Huang and Xueqin Wang(2014). Conditional distance correlation 
sure independence screening for ultrahigh dimensional data. \emph{Submitted to Biostatistics}.


Wand, M.P. & Jones, M.C. (1995) Kernel Smoothing. {Chapman & Hall}.
}
\author{
Canhong Wen, Wenliang Pan, Mian Huang, and Xueqin Wang.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{cdcor}}, \code{\link{cdcor.ada}}
}
\examples{
 set.seed(1)
 # load the distance correlation for comparison
 require(energy) 

 ## independent case
 x <- rnorm(100)
 y <- rnorm(100)
 z <- rnorm(100)
 bw(x,y,z)
 cdcor.ada(x,y,z)
 dcor(x,y)

 ## conditional dependent case
 x <- rnorm(100)
 y <- x + 0.1*rnorm(100)
 z <- rnorm(100)
 bw(x,y,z)
 cdcor.ada(x,y,z)
 dcor(x,y)

 ## conditional independent case: x and z are correlated
 require(MASS)
 data <- mvrnorm(100,rep(0,2),matrix(c(1,0.8,0.8,1),2,2))
 x <- data[,1]
 z <- data[,2]
 y <- z + 0.1*rnorm(100)
 bw(x,y,z)
 cdcor.ada(x,y,z)
 dcor(x,y)
}
\keyword{ bandwidth selection}
\keyword{ plug-in }

