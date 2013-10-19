\name{cMackSkil}
\alias{cMackSkil}
\title{
Computes a critical value for the Mack-Skillings MS distribution.
}
\description{
This function computes the critical value for the Mack-Skillings MS distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.  
}
\usage{
cMackSkil(alpha,k,n,c, method=NA, n.mc=10000)
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
  \item{c}{
  A numeric value indicating the number of replications for each treatment-block combination.
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
\item{c}{number of replications}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact")}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.9 Determination of Niacin in Bran Flakes
cMackSkil(.0501,4,3,3)
##Hollander-Wolfe-Chicken Chapter 7 Comment 72
cMackSkil(.0502,4,4,3)
}
\keyword{Mack-Skillings}