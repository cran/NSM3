\name{pNWWM}
\alias{pNWWM}
\title{
Nemenyi, Wilcoxon-Wilcox, Miller
}
\description{
Function to compute the P-value for the observed Nemenyi, Wilcoxon-Wilcox, Miller R* statistic.
}
\usage{
pNWWM(x,b=NA,trt=NA,method=NA, n.mc=10000)
}
\arguments{
  \item{x}{Either a matrix or a vector containing the data, with control assumed to be the first group.}
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

\code{pNWWM(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{pNWWM(x=c(1,2,3,4,5,6),b=c(1,1,2,2,3,3),trt=c(1,2,1,2,1,2))}
}
\value{
Returns a list with "NSM3Ch7MCp" class containing the following components:
\item{k}{number of treatments (including the control)}
\item{n}{number of blocks}
\item{obs.stat}{the observed R* statistic for each treatment vs. control comparison}
\item{p.val}{upper tail P-value corresponding to each of the k-1 observed R* statistics}
}
\author{
Grant Schneider
}
\note{
The data group containing the treatment values should be entered as the first group.
}
\examples{
##Hollander-Wolfe-Chicken Example 7.4 Stuttering Adaptation
adaptation.scores<-matrix(c(57,59,44,51,43,49,48,56,44,50,44,50,70,42,58,54,38,48,38,48,50,53,53,
56,37,58,44,50,58,48,60,58,60,38,48,56,51,56,44,44,50,54,50,40,50,50,56,46,74,57,74,48,48,44),
ncol=3,dimnames = list(1 : 18,c("No Shock", "Shock Following", "Shock During")))

#pNWWM(adaptation.scores)
pNWWM(adaptation.scores,n.mc=2500)

}
\keyword{Nemenyi}
\keyword{Wilcoxon-Wilcox}
\keyword{Miller}