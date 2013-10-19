\name{CorrUpperBound}
\alias{CorrUpperBound}
\title{
Computes the upper bound for the null correlation between two overlapping signed rank statistics.
}
\description{
This function is based on the computations in Hollander (1967).
}
\usage{
CorrUpperBound(n)
}
\arguments{
  \item{n}{number of observations}
}
\value{
Returns a numeric value indicating the upper bound.
}
\references{
Hollander, Myles. "Rank tests for randomized blocks when the alternatives have an a priori ordering." The Annals of Mathematical Statistics (1967): 867-877.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 7.12 Effect of Weight on Forearm Tremor Frequency
CorrUpperBound(6)
}
\keyword{Hollander}
\keyword{Correlation Upper Bound}