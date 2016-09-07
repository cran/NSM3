\name{cFrd}
\alias{cFrd}
\title{
Computes a critical value for the Friedman, Kendall-Babington Smith S distribution.
}
\description{
This function computes the critical value for the Friedman, Kendall-Babington Smith S distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level. The method used to compute the distribution is from the reference by Van de Wiel, Bucchianico, and Van der Laan.
}
\usage{
cFrd(alpha, k, n, method=NA, n.mc=10000, return.full.distribution=FALSE)
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
  \item{return.full.distribution}{
  If TRUE, and the method used is not asymptotic, the entire probability mass function of S will be returned.
}
}
\value{
Returns a list with "NSM3Ch7c" class containing the following components:
\item{k}{number of treatments}
\item{n}{number of blocks}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
\item{full.distribution}{probability mass function of S}
}
\references{
Van de Wiel, M. A., A. Di Bucchianico, and P. Van der Laan. "Symbolic computation and exact distributions of nonparametric test statistics." Journal of the Royal Statistical Society: Series D (The Statistician) 48.4 (1999): 507-516.
}
\author{
Grant Schneider
}
\seealso{
The \code{coin} package.
}
\examples{
##Hollander-Wolfe-Chicken Example 7.1 Rounding First Base
#cFrd(0.01,3,22,"Exact")
cFrd(0.01,3,22,n.mc=5000)
cFrd(0.01,3,22,"Asymptotic")
}
\keyword{Friedman}
\keyword{Kendall-Babington}
\keyword{Smith}