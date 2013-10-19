\name{pNDWol}
\alias{pNDWol}
\title{
Nemenyi, Damico-Wolfe
}
\description{
Function to compute the P-value for the observed Nemenyi, Damico-Wolfe Y statistic.
}
\usage{
pNDWol(x,g=NA,method=NA, n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing the data.}
 \item{g}{If x is a vector, g is a required vector of group labels. Otherwise, not used.}
   \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\value{
Returns a list with "NSM3Ch6MCp" class containing the following components:
\item{n}{number of observations in the k data groups, with the first group representing the control}
\item{obs.stat}{the observed Y statistic for each treatment vs. control comparison}
\item{p.val}{upper tail P-value corresponding to each of the k-1 observed Y statistics}
}
\author{
Grant Schneider
}
\note{
The data group containing the treatment values should be entered as the first group.
}
\examples{
##Hollander-Wolfe-Chicken Example 6.8 Motivational Effect of Knowledge of Performance
motivational.effect<-list(no.Info = c(40, 35, 38, 43, 44, 41), 
rough.Info = c(38, 40, 47, 44, 40, 42), 
accurate.Info = c(48, 40, 45, 43, 46, 44))

pNDWol(motivational.effect,method="Asymptotic")
pNDWol(motivational.effect,method="Monte Carlo")
}
\keyword{Nemenyi}
\keyword{Damico-Wolfe}