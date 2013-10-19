\name{pHayStonLSA}
\alias{pHayStonLSA}
\title{
Hayter-Sone LSA
}
\description{
Function to compute the upper tail probability of the Hayter-Stone W asymptotic distribution for a given cutoff.
}
\usage{
pHayStonLSA(h,k,delta=.001)
}
\arguments{
  \item{h}{Cutoff used to calculate the P-value.}
  \item{k}{Number of groups.}
  \item{delta}{Defines the fineness of the grid used to calculate the asymptotic distribution of W.}
}
\value{
Returns the asymptotic upper tail P-value.
}
\author{
Grant Schneider
}
\examples{
pHayStonLSA(2.491,3)
pHayStonLSA(4.112,4)
}
\keyword{Hayter-Stone}
\keyword{LSA}