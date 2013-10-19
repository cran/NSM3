\name{dmrl.mc}
\alias{dmrl.mc}
\title{
Hollander-Proschan
}
\description{
Function to compute the Monte Carlo or asymptotic P-value for the observed Hollander-Proschan V' statistic.
}
\usage{
dmrl.mc(x, alternative = "two.sided", exact=FALSE,
        min.reps = 100, max.reps = 1000, delta = 10^-3)
}
\arguments{
   \item{x}{a vector of data of length n}
  \item{alternative}{the direction of the alternative hypothesis.  The choices are two.sided, dmrl, and imrl with the default value being two.sided.}
  \item{exact}{TRUE/FALSE value that determines whether the exact test or the large sample approximation is used if n >= 9.  If n < 9 the exact test is used.    The default value is FALSE, so the large sample approximation will be used unless specified not to.  This is the same large sample approximation as epstein()}
  \item{min.reps}{the minimum number of repetitions for the Monte Carlo Approximation}
\item{max.reps}{the maximum number of reps for the Monte Carlo Approximation. If the maximum number of reps has been reached, and the probability has not converged, a warning is given.}    
\item{delta}{the measure of accuracy for the convergence.  If the probability converges to within delta, the Monte Carlo procedure stops before reaching the maximum number of reps.}
}
\value{
The function returns a list with two elements:
\item{V}{the value of the dmrl statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
ex11.1<-c(42, 43, 51, 61, 66, 69, 71, 81, 82, 82)
dmrl.mc(ex11.1, alt="dmrl", exact=TRUE)

}
\keyword{Hollander-Proschan}
\keyword{DMRL}
\keyword{Monte Carlo}