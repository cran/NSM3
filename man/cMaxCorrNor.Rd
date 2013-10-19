\name{cMaxCorrNor}
\alias{cMaxCorrNor}
\title{
Quantile function for the maximum of k N(0,1) random variables with common correlation rho.
}
\description{
Uses the integrate function based on the method proposed in Gupta, Panchapakesan and Sohn (1983).
}
\usage{
cMaxCorrNor(alpha,k,rho)
}
\arguments{
  \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{k}{Number of random variables.}
  \item{rho}{Common correlation between the random variables.}
}
\value{
Returns the upper tail cutoff at or immediately below the user-specified alpha.
}
\references{
Gupta, Shanti S., S. Panchapakesan, and Joong K. Sohn. "On the distribution of the studentized maximum of equally correlated normal random variables." Communications in Statistics-Simulation and Computation 14.1 (1985): 103-135.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Section 7.4 LSA
cMaxCorrNor(.04584,4,.5)
##Hollander-Wolfe-Chicken Section 7.14
cMaxCorrNor(.02337,5,.3)
##Hollander-Wolfe-Chicken Example 7.14
cMaxCorrNor(.10,5,.452)
}
\keyword{Maximum Correlated Normal}