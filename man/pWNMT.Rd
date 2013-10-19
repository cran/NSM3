\name{pWNMT}
\alias{pWNMT}
\title{
Wilcoxon, Nemenyi, McDonald-Thompson
}
\description{
Function to compute the P-value for the observed Wilcoxon, Nemenyi, McDonald-Thompson R statistic.
}
\usage{
pWNMT(x,b=NA,trt=NA,method=NA, n.mc=10000)
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
\code{pWNMT(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{pWNMT(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7MCp" class containing the following components:
\item{k}{number of treatments}
\item{n}{number of blocks}
\item{obs.stat}{the observed R* statistic for each of the k*(k-1)/2 comparisons}
\item{p.val}{upper tail P-value corresponding to each observed R statistic}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.3 Rounding First Base
RoundingTimes<-matrix(c(5.40, 5.50, 5.55, 5.85, 5.70, 5.75, 5.20, 5.60, 5.50, 5.55, 5.50, 5.40,
5.90, 5.85, 5.70, 5.45, 5.55, 5.60, 5.40, 5.40, 5.35, 5.45, 5.50, 5.35, 5.25, 5.15, 5.00, 5.85,
5.80, 5.70, 5.25, 5.20, 5.10, 5.65, 5.55, 5.45, 5.60, 5.35, 5.45, 5.05, 5.00, 4.95, 5.50, 5.50,
5.40, 5.45, 5.55, 5.50, 5.55, 5.55, 5.35, 5.45, 5.50, 5.55, 5.50, 5.45, 5.25, 5.65, 5.60, 5.40,
5.70, 5.65, 5.55, 6.30, 6.30, 6.25),nrow = 22,byrow = TRUE,dimnames = list(1 : 22,
c("Round Out", "Narrow Angle", "Wide Angle")))

pWNMT(RoundingTimes,n.mc=2500)
}
\keyword{Wilcoxon}
\keyword{Nemenyi}
\keyword{McDonald-Thompson}