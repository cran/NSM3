\name{klefsjo.ifra.mc}
\alias{klefsjo.ifra.mc}
\title{
Function to compute the Monte Carlo P-value for the observed Klefsjo's B* statistic.
}
\description{
This is the Monte Carlo approximation to the function "klefsjo.ifra".
}
\usage{
klefsjo.ifra.mc(x, alternative = "two.sided", exact=FALSE,
                min.reps = 100, max.reps = 1000, delta = 10^-3)
}
\arguments{
   \item{x}{a vector of data of length n}
  \item{alternative}{the direction of the alternative hypothesis.  The choices are two.sided, ifra and dfra with the default value being two.sided.}
  \item{exact}{TRUE/FALSE value that determines whether the exact test or the large sample approximation is used if n >= 9.  If n < 9 the exact test is used.    The default value is FALSE, so the large sample approximation will be used unless specified not to.  This is the same large sample approximation as epstein()}
  \item{min.reps}{the minimum number of repetitions for the Monte Carlo Approximation}
\item{max.reps}{the maximum number of reps for the Monte Carlo Approximation. If the maximum number of reps has been reached, and the probability has not converged, a warning is given.}    
\item{delta}{the measure of accuracy for the convergence.  If the probability converges to within delta, the Monte Carlo procedure stops before reaching the maximum number of reps.}
}
\value{
The function returns a list with two elements:
\item{B.star}{the value of the Klefsjo statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
temp.data<-c(0.33925023, 0.84005767, 0.29066189, 1.95163010, 0.74536608, 0.16714902, 0.06950791, 
1.14919291,1.93210982, 1.06006126, 0.14651009, 0.28776282, 0.72242750, 1.02227211, 1.71243334)
klefsjo.ifra.mc(temp.data, exact=TRUE)
}
\keyword{Klefsjo}
\keyword{IFRA}
\keyword{Monte Carlo}