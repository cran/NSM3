\name{pRangeNor}
\alias{pRangeNor}
\title{
Function to compute the upper-tail probability of the range of k independent N(0,1) random variables for a given cutoff.
}
\description{
Uses the integrate function based on the method proposed in Harter (1960).
}
\usage{
pRangeNor(x,k)
}
\arguments{
 \item{x}{Cutoff at which the upper-tail P-value is to be calculated.}
  \item{k}{Number of independent Normal random variables.}
}
\value{
Returns the upper tail probability at the user-specified cutoff.
}
\references{
Harter, H. Leon. "Tables of range and studentized range." The Annals of Mathematical Statistics (1960): 1122-1147.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.3 Rounding First Base
pRangeNor(4.121,3)

##Hollander-Wolfe-Chicken Example 7.7 Chemical Toxicity
pRangeNor(4.171,7)
}
\keyword{Range of Independent Normal}