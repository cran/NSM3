\name{nb.mc}
\alias{nb.mc}
\title{
Function to compute the Monte Carlo P-value for the observed Hollander-Proschan T statistic.
}
\description{
This is the Monte Carlo approximation to the newbet function.
}
\usage{
nb.mc(x, alternative = "two.sided", exact=FALSE, 
      min.reps = 100, max.reps = 1000, delta = 10^-3)
}
\arguments{
   \item{x}{a vector of data of length n}
  \item{alternative}{the direction of the alternative hypothesis.  The choices are two.sided, nbu, and nwu with the default value being two.sided.}
  \item{exact}{TRUE/FALSE value that determines whether the exact test or the large sample approximation is used if n >= 9.  If n < 9 the exact test is used.    The default value is FALSE, so the large sample approximation will be used unless specified not to.  This is the same large sample approximation as epstein()}
  \item{min.reps}{the minimum number of repetitions for the Monte Carlo Approximation}
\item{max.reps}{the maximum number of reps for the Monte Carlo Approximation. If the maximum number of reps has been reached, and the probability has not converged, a warning is given.}    
\item{delta}{the measure of accuracy for the convergence.  If the probability converges to within delta, the Monte Carlo procedure stops before reaching the maximum number of reps.}

}
\value{
The function returns a list with two elements:
\item{T}{the value of the Hollander-Proschan statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
table11.4<-c(194,15,41,29,33,181)
nb.mc(table11.4, alt="nbu")
}
\keyword{Hollander-Proschan}
\keyword{Monte Carlo}