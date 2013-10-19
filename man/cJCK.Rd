\name{cJCK}
\alias{cJCK}
\title{
Computes a critical value for the Jonckheere-Terpstra J distribution.
}
\description{
This function computes the critical value for the Jonckheere-Terpstra J distribution at (or typically in the "Exact" case, close to) the given alpha level. The function takes advantage of Harding's (1984) algorithm to quickly generate the distribution.
}
\usage{
cJCK(alpha, n, method=NA, n.mc=10000)
}
\arguments{
  \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{n}{
  A vector of numeric values indicating the size of each of the k data groups.
}
  \item{method}{
  Either "Exact" or "Asymptotic", indicating the desired distribution. When method=NA, if sum(n)<=200, the "Exact" method will be used to compute the J distribution. Otherwise, the "Asymptotic" method will be used.
}
  \item{n.mc}{
  Not used.  Only included for standardization with other critical value procedures in the NSM3 package.
}
}
\value{
Returns a list with "NSM3Ch6c" class containing the following components:
\item{n}{number of observations in the k data groups}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact")}
}
\references{
Harding, E. F. "An efficient, minimal-storage procedure for calculating the Mann-Whitney U, generalized U and similar distributions." Applied statistics (1984): 1-6.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 6.2 Motivational Effect of Knowledge of Performance
cJCK(.0490, c(6,6,6),"Exact")
cJCK(.0490, c(6,6,6),"Monte Carlo")
cJCK(.0231, c(6,6,6),"Exact")
}
\keyword{Jonckheere-Terpstra}