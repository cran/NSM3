\name{cSDCFlig}
\alias{cSDCFlig}
\title{
Computes a critical value for the Dwass, Steel, Critchlow-Fligner W distribution.
}
\description{
This function computes the critical value for the Dwass, Steel, Critchlow-Fligner W distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cSDCFlig(alpha, n, method=NA, n.mc=10000)
}
\arguments{
   \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{n}{
  A vector of numeric values indicating the size of each of the k data groups.
}
  \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\value{
Returns a list with "NSM3Ch6c" class containing the following components:
\item{n}{number of observations in the k data groups}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Chapter 6 Comment 55
#cSDCFlig(.0331, c(3, 5, 7),n.mc=10000)
cSDCFlig(.0331, c(3, 5, 7),n.mc=2500)

##Another example
#cSDCFlig(alpha=0.05,n=rep(4,3),method="Exact")
cSDCFlig(alpha=0.05,n=rep(4,3),method="Monte Carlo",n.mc=2500)
#cSDCFlig(alpha=0.05,n=rep(4,3),method="Asymptotic")
}
\keyword{Dwass}
\keyword{Steel}
\keyword{Critchlow-Fligner}