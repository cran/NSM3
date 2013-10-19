\name{pJCK}
\alias{pJCK}
\title{
Function to compute the P-value for the observed Jonckheere-Terpstra J statistic.
}
\description{
This function computes the observed J statistic for the given data and corresponding P-value. When there are no ties in the data, the function takes advantage of Harding's (1984) algorithm to quickly generate the exact distribution of J.
}
\usage{
pJCK(x,g=NA,method=NA, n.mc=10000)
}
\arguments{

 \item{x}{Either a list or a vector containing the data.}
 \item{g}{If x is a vector, g is a required vector of group labels. Otherwise, not used.}
  \item{method}{Either "Exact", "Monte Carlo", or "Asymptotic", indicating the desired distribution. When method=NA and ties are not present, "Exact" will be used. When method=NA and ties are present, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
 \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the groups of data can be entered in either of two ways. For data a=1,2 and b=3,4,5 the following are equivalent:

\code{pJCK(x=list(c(1,2),c(3,4,5)))}
\code{pJCK(x=c(1,2,3,4,5),g=c(1,1,2,2,2))}

}
\value{
Returns a list with "NSM3Ch6p" class containing the following components:
\item{n}{a vector containing the number of observations in each of the data groups}
\item{obs.stat}{the observed J statistic}
\item{p.val}{upper tail P-value}
}
\references{
Harding, E. F. "An efficient, minimal-storage procedure for calculating the Mann-Whitney U, generalized U and similar distributions." Applied statistics (1984): 1-6.
}
\author{
Grant Schneider
}

\examples{
##Hollander-Wolfe-Chicken Example 6.2 Motivational Effect of Knowledge of Performance
motivational.effect<-list(no.Info=c(40,35,38,43,44,41),rough.Info=c(38,40,47,44,40,42),
                          accurate.Info=c(48,40,45,43,46,44))
#pJCK(motivational.effect,method="Monte Carlo")
pJCK(motivational.effect,method="Asymptotic")
}
\keyword{Jonckheere-Terpstra}