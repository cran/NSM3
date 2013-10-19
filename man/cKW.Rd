\name{cKW}
\alias{cKW}
\title{
Computes a critical value for the Kruskal-Wallis H distribution.
}
\description{
This function computes the critical value for the Kruskal-Wallis H distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cKW(alpha,n, method=NA, n.mc=10000)
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
##Hollander-Wolfe-Chicken Example 6.1 Half-Time of Mucociliary Clearance
#cKW(0.0503,c(5,4,5),"Exact")
cKW(0.7147,c(5,4,5),"Asymptotic")
cKW(0.7147,c(5,4,5),"Monte Carlo",n.mc=20000)
}
\keyword{Kruskal-Wallis}