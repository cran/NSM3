\name{e.mc}
\alias{e.mc}
\title{
Function to compute the Monte Carlo P-value for the observed Epstein E statistic
}
\description{
This is the Monte Carlo approximation to the function "epstein".
}
\usage{
e.mc(x, alternative = "two.sided", exact=FALSE,
     min.reps = 100, max.reps = 1000, delta = 10^-3)
}
\arguments{
  \item{x}{a vector of data of length n}
  \item{alternative}{the direction of the alternative hypothesis.  The choices are two.sided, ifr and dfr with the default value being two.sided.}
  \item{exact}{TRUE/FALSE value that determines whether the exact test or the large sample approximation is used if n >= 9.  If n < 9 the exact test is used.    The default value is FALSE, so the large sample approximation will be used unless specified not to.  This is the same large sample approximation as epstein()}
  \item{min.reps}{the minimum number of repetitions for the Monte Carlo Approximation}
\item{max.reps}{the maximum number of reps for the Monte Carlo Approximation. If the maximum number of reps has been reached, and the probability has not converged, a warning is given.}    
\item{delta}{the measure of accuracy for the convergence.  If the probability converges to within delta, the Monte Carlo procedure stops before reaching the maximum number of reps.}  
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
Ep <- e.mc(ex11.1, alt="ifr", exact=TRUE)
Ep$E
Ep$p

#Large Sample Approximation
Ep.lsa <- e.mc(ex11.1, alt="ifr")

table11.2<-c(487, 18, 100, 7, 98, 5, 85, 91, 43, 230, 3, 130)
Ep=e.mc(table11.2,alt="i", exact=TRUE)
#Failing to converge
Ep=e.mc(table11.2,alt="i", exact=TRUE, min.reps=5, max.reps=5)
}

\keyword{Epstein}
\keyword{Monte Carlo}
