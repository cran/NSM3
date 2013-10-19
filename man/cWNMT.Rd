\name{cWNMT}
\alias{cWNMT}
\title{
Computes a critical value for the Wilcoxon, Nemenyi, McDonald-Thompson R distribution.
}
\description{
This function computes the critical value for the Wilcoxon, Nemenyi, McDonald-Thompson R distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cWNMT(alpha, k, n, method=NA, n.mc=10000)
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
##Hollander-Wolfe-Chicken Example 7.3 Rounding First Base
#cWNMT(.047, 3, 15)
cWNMT(.047, 3, 15,n.mc=5000)

##Chapter 7 Comment 26
#cWNMT(.083, 4, 2)
cWNMT(.083, 4, 2,n.mc=5000)
}
\keyword{Wilcoxon}
\keyword{Nemenyi}
\keyword{McDonald-Thompson}