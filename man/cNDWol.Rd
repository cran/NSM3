\name{cNDWol}
\alias{cNDWol}
\title{
Function to compute a critical value for the Nemenyi, Damico-Wolfe Y distribution.
}
\description{
This function computes the critical value for the Nemenyi, Damico-Wolfe Y distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cNDWol(alpha,n, method=NA, n.mc=10000)
}
\arguments{
\item{alpha}{
  A numeric value between 0 and 1.
}
  \item{n}{
  A vector of numeric values indicating the size of each of the k data groups, with the first element indicating the treatment group size.
}
   \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\value{
Returns a list with "NSM3Ch6MCc" class containing the following components:
\item{n}{number of observations in the k data groups}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 6.8 Motivational Effect of Knowledge of Performance
cNDWol(.0554, c(6, 6, 6),"Monte Carlo")
cNDWol(.0554, c(6, 6, 6),"Monte Carlo",n.mc=25000)
cNDWol(.0371, c(6, 6, 6),"Monte Carlo")
}
\keyword{Nemenyi}
\keyword{Damico-Wolfe}