\name{cSkilMack}
\alias{cSkilMack}
\title{
Computes a critical value for the Skillings-Mack SM distribution.
}
\description{
This function computes the critical value for the Skillings-Mack SM distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.  
}
\usage{
cSkilMack(alpha, obs.mat, method = NA, n.mc = 10000)
}
\arguments{
 \item{alpha}{A numeric value between 0 and 1.}
  \item{obs.mat}{The incidence matrix, explained below.}
  \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The incidence matrix, obs.mat, will be an n x k matrix of ones and zeroes, which indicate where the data are observed and unobserved, respectively.}
\value{
Returns a list with "NSM3Ch7c" class containing the following components:
 \item{k}{number of treatments}
 \item{n}{number of blocks}
 \item{s}{number of treatments per block}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
}
\author{
Grant Schneider
}
\note{
The syntax of this procedure differs from the others in the NSM3 package due to the fact that the distribution is calculated conditionally on the pattern of missingness. We therefore require obs.mat, the incidence matrix.
}

\examples{
##Hollander, Wolfe, Chicken Example 7.8 Effect of Rhythmicity of a Metronome on Speech Fluency
obs.mat<-matrix(c(rep(1,10),0,rep(1,13)),ncol=3,byrow=TRUE)
#cSkilMack(.01,obs.mat)
cSkilMack(.01,obs.mat,n.mc=5000)
}
\keyword{Skillings-Mack}