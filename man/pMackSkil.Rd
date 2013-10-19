\name{pMackSkil}
\alias{pMackSkil}
\title{
Mack-Skillings
}
\description{
Function to compute the P-value for the observed Mack-Skillings MS statistic.
}
\usage{
pMackSkil(x,b=NA,trt=NA,method=NA,n.mc=10000)
}
\arguments{
  \item{x}{Either a 3 dimensional array or a vector containing the data.}
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

\code{pMackSkil(x=array(c(1,2,3,4,5,6),dim=c(1,2,3))}
\code{pMackSkil(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7p" class containing the following components:
\item{k}{number of treatments in the data}
\item{n}{number of blocks in the data}
\item{c}{number of repetitions for each treatment and block combination}
\item{obs.stat}{the observed MS statistic}
\item{p.val}{upper tail P-value}
}
\author{
Grant Schneider
}
\examples{
##Hollander, Wolfe, Chicken Example 7.9 Determination of Niacin in Bran Flakes
niacin<-array(dim=c(3,4,3))
niacin[,,1]<-c(7.58,7.87,7.71,8,8.27,8,7.6,7.3,7.82,8.03,7.35,7.66)
niacin[,,2]<-c(11.63,11.87,11.4,12.2,11.7,11.8,11.04,11.5,11.49,11.5,10.10,11.7)
niacin[,,3]<-c(15,15.92,15.58,16.6,16.4,15.9,15.87,15.91,16.28,15.1,14.8,15.7)
}
\keyword{Mack-Skillings}