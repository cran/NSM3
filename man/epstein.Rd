\name{epstein}
\alias{epstein}
\title{
Epstein
}
\description{
Function to compute the P-value for the observed Epstein E statistic
}
\usage{
epstein(x, alternative = "two.sided", exact=FALSE)
}
\arguments{
  \item{x}{a vector of data of length n}
  \item{alternative}{the direction of the alternative hypothesis.  The choices are two.sided, ifr (for increasing failure rate) and dfr (for decreasing failure rate) with the default value being two.sided.}
  \item{exact}{TRUE/FALSE value that determines whether the exact test or the large sample approximation is used if n >= 9. If n < 9 the exact test is used. The default value is FALSE, so the large sample approximation will be used unless specified not to.}    

}
\value{
The function returns a list with two elements:
\item{E}{the value of the Epstein statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
ex11.1<-c(42, 43, 51, 61, 66, 69, 71, 81, 82, 82)
Ep <- epstein(ex11.1, alt="ifr", exact=TRUE)
Ep$E
Ep$p

#Large Sample Approximation
Ep.lsa <- epstein(ex11.1, alt="ifr")

}
\keyword{epstein}