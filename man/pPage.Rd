\name{pPage}
\alias{pPage}
\title{
Page
}
\description{
Function to compute the P-value for the observed Page L statistic.
}
\usage{
pPage(x,b=NA,trt=NA,method=NA, n.mc=10000)
}
\arguments{
  \item{x}{Either a matrix or a vector containing the data.}
  \item{b}{If x is a vector, b is a required vector of block labels. Otherwise, not used.}
  \item{trt}{If x is a vector, trt is a required vector of treatment labels. Otherwise, not used.}
  \item{method}{Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the data can be entered in either of two ways. The following are equivalent:

\code{pPage(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{pPage(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7p" class containing the following components:
\item{k}{number of treatments in the data}
\item{n}{number of blocks in the data}
\item{obs.stat}{the observed L statistic}
\item{p.val}{upper tail P-value}
}

\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.2 Breaking Strength of Cotton Fibers
strength.index<-matrix(c(7.46, 7.68, 7.21, 7.17, 7.57, 7.80, 7.76, 7.73, 7.74, 8.14, 8.15,
7.87, 7.63, 8.00, 7.93),byrow=FALSE,ncol=5)

#pPage(strength.index,method="Exact")
pPage(strength.index,method="Monte Carlo")
}
\keyword{Page}