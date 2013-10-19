\name{pSkilMack}
\alias{pSkilMack}
\title{
Skillings-Mack
}
\description{
Function to compute the P-value for the observed Skillings-Mack SM statistic.
}
\usage{
pSkilMack(x, b = NA, trt = NA, method = NA, n.mc = 10000)
}
\arguments{
  \item{x}{Either a matrix or a vector containing the data.}
  \item{b}{If x is a vector, b is a required vector of block labels. Otherwise, not used.}
  \item{trt}{If x is a vector, trt is a required vector of treatment labels. Otherwise, not used.}
  \item{method}{Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}

\details{
The data entry is intended to be flexible, so that the data can be entered in either of two ways. The following are equivalent:
\code{pSkilMack(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{pSkilMack(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7p" class containing the following components:
\item{k}{number of treatments in the data}
\item{n}{number of blocks in the data}
\item{s}{number of treatments per block}
\item{obs.stat}{the observed D statistic}
\item{p.val}{upper tail P-value}
}
\author{
Grant Schneider
}
\examples{
##Hollander, Wolfe, Chicken Example 7.8 Effect of Rhythmicity of a Metronome on Speech Fluency
rhythmicity<-matrix(c(3, 5, 15, 1, 3, 18, 5, 4, 21, 2, NA, 6, 0, 2, 17, 0, 2, 10, 0, 3, 8,
0, 2, 13),ncol=3,byrow=TRUE)
#pSkilMack(rhythmicity)
pSkilMack(rhythmicity,n.mc=5000)
}
\keyword{Skillings-Mack}