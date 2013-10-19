\name{klefsjo.ifr}
\alias{klefsjo.ifr}
\title{
Klefsjo's IFR
}
\description{
Function to compute the P-value for the observed Klefsjo's A* statistic.
}
\usage{
klefsjo.ifr (x, alternative = "two.sided", exact=FALSE)
}
\arguments{
\item{x}{a vector of data of length n}
\item{alternative}{the direction of the alternative hypothesis.  The choices are two.sided, ifr and dfr with the default value being two.sided.}  
\item{exact}{TRUE/FALSE value that determines whether the exact test or the large sample approximation is used if n >= 9.  If n < 9 the exact test is used. The default value is FALSE, so the large sample approximation will be used unless specified not to.}
}
\details{
If the sample size is too large to allow for an exact value, due to duplicate coefficients, a note will be displayed and the large sample approximation will be used. 
}
\value{
The function returns a list with two elements:
\item{A.star}{the value of the Klefsjo statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
velocity<-c(12.8, 12.9, 13.3, 13.4, 13.7, 13.8, 14.5)
klefsjo.ifr(velocity)

#Example of forced Large Sample Approximation
tb<-c(43, 45, 53, 56, 56, 57, 58, 66, 67, 73, 74, 79, 80, 80, 81, 81, 81, 82, 83, 83, 84, 88,
89,  91,  91,  92,  92,  97,  99,  99, 100, 100, 101, 102, 102, 102, 103, 104, 107, 108, 109,
113, 114, 118, 121, 123, 126, 128, 137, 138, 139, 144, 145, 147, 156, 162, 174, 178, 179, 184,
191, 198, 211, 214, 243, 249, 329, 380, 403, 511, 522, 598)
klefsjo.ifr(tb, exact=TRUE)
}
\keyword{Klefsjo}
\keyword{IFR}