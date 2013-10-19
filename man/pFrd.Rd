\name{pFrd}
\alias{pFrd}
\title{
Function to compute the P-value for the observed Friedman, Kendall-Babington Smith S statistic.
}
\description{
The method used to compute the P-value is from the reference by Van de Wiel, Bucchianico, and Van der Laan.
}
\usage{
pFrd(x,b=NA,trt=NA,method=NA, n.mc=10000)
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

\code{pFrd(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{pFrd(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7p" class containing the following components:
\item{k}{number of treatments in the data}
\item{n}{number of blocks in the data}
\item{obs.stat}{the observed D statistic}
\item{p.val}{upper tail P-value}
}

\references{
Van de Wiel, M. A., A. Di Bucchianico, and P. Van der Laan. "Symbolic computation and exact distributions of nonparametric test statistics." Journal of the Royal Statistical Society: Series D (The Statistician) 48.4 (1999): 507-516.
}
\author{
Grant Schneider
}
\seealso{
Also see the \code{coin} package.
}
\examples{
##Hollander-Wolfe-Chicken Example 7.1 Rounding First Base
rounding.times<-matrix(c(5.40, 5.50, 5.55,
                         5.85, 5.70, 5.75,
                         5.20, 5.60, 5.50,
                         5.55, 5.50, 5.40,
                         5.90, 5.85, 5.70,
                         5.45, 5.55, 5.60,
                         5.40, 5.40, 5.35,
                         5.45, 5.50, 5.35,
                         5.25, 5.15, 5.00,
                         5.85, 5.80, 5.70,
                         5.25, 5.20, 5.10,
                         5.65, 5.55, 5.45,
                         5.60, 5.35, 5.45,
                         5.05, 5.00, 4.95,
                         5.50, 5.50, 5.40,
                         5.45, 5.55, 5.50,
                         5.55, 5.55, 5.35,
                         5.45, 5.50, 5.55,
                         5.50, 5.45, 5.25,
                         5.65, 5.60, 5.40,
                         5.70, 5.65, 5.55,
                         6.30, 6.30, 6.25),ncol=3,byrow=TRUE)
#pFrd(rounding.times,n.mc=20000)
pFrd(rounding.times,n.mc=2000)
}
\keyword{Friedman}
\keyword{Kendall-Babington Smith}