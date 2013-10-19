\name{cHaySton}
\alias{cHaySton}
\title{
Computes a critical value for the Hayter-Stone W* distribution.
}
\description{
This function computes the critical value for the Hayter-Stone W* distriburion at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.
}
\usage{
cHaySton(alpha,n, method=NA, n.mc=10000)
}
\arguments{
   \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{n}{
  A vector (of length 2 or greater) indicating the sizes of the data groups.
}
  \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the         desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The Asymptotic distribution requires that all group sizes are equal. If method="Asymptotic" and there are different group sizes in n, method="Monte Carlo" will be used.
}
\value{
Returns a list with "NSM3Ch6MCc" class containing the following components:
\item{n}{data group sizes}
\item{num.comp}{number of multiple comparisons to be made (based on the length of n)}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")
}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 6.7 Motivational Effect of Knowledge of Performance:
#cHaySton(.0553,rep(6,3),"Monte Carlo")
cHaySton(.05,c(6,6,6),"Asymptotic")
}
\keyword{Hayter-Stone}