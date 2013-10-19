\name{cRangeNor}
\alias{cRangeNor}
\title{
Quantile function for the range of k independent N(0,1) random variables.
}
\description{
Uses the integrate function based on the method proposed in Harter (1960).
}
\usage{
cRangeNor(alpha,k)
}
\arguments{
   \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{k}{Number of independent Normal random variables.}
}
\value{
Returns the upper tail cutoff at or immediately below the user-specified alpha.
}
\references{
Harter, H. Leon. "Tables of range and studentized range." The Annals of Mathematical Statistics (1960): 1122-1147.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.3 Rounding First Base
cRangeNor(.01, 3)

##Hollander-Wolfe-Chicken Example 7.7 Chemical Toxicity
cRangeNor(.05, 7)
}
\keyword{Range of Independent Normal}