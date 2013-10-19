\name{cPage}
\alias{cPage}
\title{
Function to compute a critical value for the Page L distribution.
}
\description{
This function computes the critical value for the Page L distriburion at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cPage(alpha, k, n, method=NA, n.mc=10000)
}
\arguments{
  \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{k}{
  A numeric value indicating the number of treatments.
}
  \item{n}{
  A numeric value indicating the number of blocks.
}
  \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\value{
Returns a list with "NSM3Ch7c" class containing the following components:
\item{k}{number of treatments}
\item{n}{number of blocks}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.2 Breaking Strength of Cotton Fibers
#cPage(.0097, 5, 3,"Exact")
cPage(.0097, 5, 3,"Monte Carlo")
}
\keyword{Page}