\name{pMaxCorrNor}
\alias{pMaxCorrNor}
\title{
Function to compute the upper tail probability of the maximum of k N(0,1) random variables with common correlation for a given cutoff.
}
\description{
Uses the integrate function based on the method proposed in Gupta, Panchapakesan and Sohn (1983).
}
\usage{
pMaxCorrNor(x,k,rho)
}
\arguments{
  \item{x}{Cutoff at which the upper-tail P-value is to be calculated.}
  \item{k}{Number of random variables.}
  \item{rho}{Common correlation between the random variables.}

}
\value{
Returns the upper tail probability at the user-specified cutoff.
}
\references{
Gupta, Shanti S., S. Panchapakesan, and Joong K. Sohn. "On the distribution of the studentized maximum of equally correlated normal random variables." Communications in Statistics-Simulation and Computation 14.1 (1985): 103-135.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Section 7.14
pMaxCorrNor(2.575,5,.3)

##Hollander-Wolfe-Chicken Example 7.14 Effect of Weight on Forearm Tremor Frequency
pMaxCorrNor(1.93,5,.452)

}
\keyword{Maximum Correlated Normal}