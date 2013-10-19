\name{cLepage}
\alias{cLepage}
\title{
Computes a critical value for the Lepage D distribution.
}
\description{
This function computes the critical value for the Lepage D distriburion at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cLepage(alpha, m, n, method=NA, n.mc=10000)
}
\arguments{
  \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{m}{
  A numeric value indicating the size of the first data group (X).
}
  \item{n}{
  A numeric value indicating the size of the second data group (Y).
}
  \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\value{
Returns a list with "NSM3Ch5c" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 5.3 Platelet Counts of Newborn Infants
cLepage(0.02,10,6,"Exact")
cLepage(0.02,10,6,"Monte Carlo")
cLepage(0.02,10,6,"Asymptotic")
}
\keyword{Lepage}