\name{cDurSkiMa}
\alias{cDurSkiMa}
\title{
Computes a critical value for the Durbin, Skillings-Mack D distribution.
}
\description{
This function computes the critical value for the Durbin, Skillings-Mack D distribution at (or typically in the "Exact" and "Monte Carlo" cases, close to) the given alpha level.  
}
\usage{
cDurSkiMa(alpha,obs.mat, method=NA, n.mc=10000)
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
The incidence matrix, obs.mat, will be an n x k matrix of ones and zeroes, which indicate where the data are observed and unobserved, respectively. Methods for finding the incidence matrix for various BIBD designs are given in the literature. While the incidence matrix will not be unique for a given (k, n, s, lambda, p) combination, the distribution of D under H0 will be the same.
}
\value{
Returns a list with "NSM3Ch7c" class containing the following components:
 \item{k}{number of treatments}
 \item{n}{number of blocks}
 \item{s}{number of treatments per block}
 \item{pp}{number of observations per treatment}
 \item{lambda}{number of times each pair of treatments occurs together within a block}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact" or "Monte Carlo")}
}
\author{
Grant Schneider
}
\note{
The syntax of this procedure differs from the others in the NSM3 package due to the fact that creating a BIBD for a given k,n,s,p,lambda is not trivial. We therefore require obs.mat, the incidence matrix.
}
\examples{
##Hollander, Wolfe, Chicken Chapter 7, comment 49 
obs.mat<-matrix(c(1,1,0,1,0,1,0,1,1),ncol=3,byrow=TRUE)
cDurSkiMa(.75,obs.mat)
}
\keyword{Durbin}
\keyword{Skillings-Mack}