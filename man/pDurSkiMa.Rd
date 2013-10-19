\name{pDurSkiMa}
\alias{pDurSkiMa}
\title{
Durbin, Skillings-Mack
}
\description{
Function to compute the P-value for the observed Durbin, Skillings-Mack D statistic.
}
\usage{
pDurSkiMa(x,b=NA,trt=NA,method=NA,n.mc=10000)
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
\code{pDurSkiMa(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{pDurSkiMa(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7p" class containing the following components:
\item{k}{number of treatments in the data}
\item{n}{number of blocks in the data}
\item{ss}{number of treatments per block}
 \item{pp}{number of observations per treatment}
 \item{lambda}{number of times each pair of treatments occurs together within a block}
\item{obs.stat}{the observed D statistic}
\item{p.val}{upper tail P-value}
}
\author{
Grant Schneider
}
\examples{
##Hollander, Wolfe, Chicken Example 7.6 Chemical Toxicity
table7.12<-matrix(nrow=7,ncol=7)
table7.12[1,c(1,2,4)]<-c(0.465,0.343,0.396)
table7.12[2,c(1,3,5)]<-c(0.602,0.873,0.634)
table7.12[3,c(3,4,7)]<-c(0.875,0.325,0.330)
table7.12[4,c(1,6,7)]<-c(0.423,0.987,0.426)
table7.12[5,c(2,3,6)]<-c(0.652,1.142,0.989)
table7.12[6,c(2,5,7)]<-c(0.536,0.409,0.309)
table7.12[7,c(4,5,6)]<-c(0.609,0.417,0.931)

pDurSkiMa(table7.12)

##or, equivalently:
x<-c(.465,.602,.423,.343,.652,.536,.873,.875,1.142,.396,.325,.609,.634,.409,.417,.987,.989,
.931,.330,.426,.309)
b<-c(1,2,4,1,5,6,2,3,5,1,3,7,2,6,7,4,5,7,3,4,6)
trt<-c(rep("A",3),rep("B",3),rep("C",3),rep("D",3),rep("E",3),rep("F",3),rep("g",3))

pDurSkiMa(x,b,trt)

}
\keyword{Durbin}
\keyword{Skillings-Mack}