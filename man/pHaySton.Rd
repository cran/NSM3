\name{pHaySton}
\alias{pHaySton}
\title{
Hayter-Stone
}
\description{
Function to compute the P-value for the observed Hayter-Stone W statistic.
}
\usage{
pHaySton(x,g=NA,method=NA,n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing the data.}
 \item{g}{If x is a vector, g is a required vector of group labels. Otherwise, not used.}
 \item{method}{Either "Exact", "Monte Carlo", or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the groups of data can be entered in either of two ways. For data a=1,2 and b=3,4,5 the following are equivalent:

\code{pHaySton(x=list(c(1,2),c(3,4,5)))}
\code{pHaySton(x=c(1,2,3,4,5),g=c(1,1,2,2,2))}

}
\value{
Returns a list with "NSM3Ch6MCp" class containing the following components:
\item{n}{a vector containing the number of observations in each of the data groups}
\item{obs.stat}{the observed W statistic for each of the k*(k-1)/2 comparisons}
\item{p.val}{upper tail P-value corresponding to each W statistic}
}
\author{
Grant Schneider
}
\examples{
##Hollander, Wolfe, Chicken Example 6.7 Motivational Effect of Knowledge of Performance:
motivational.effect<-list(no.Info = c(40, 35, 38, 43, 44, 41), rough.Info = c(38, 
40, 47, 44, 40, 42), accurate.Info = c(48, 40, 45, 43, 46, 44
))

#pHaySton(motivational.effect,method="Monte Carlo")
pHaySton(motivational.effect,method="Asymptotic")
#pHaySton(rnorm(10),rep(1:3,c(3,3,4)),method="Asymptotic")
}
\keyword{Hayter-Stone}