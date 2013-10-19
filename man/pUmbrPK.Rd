\name{pUmbrPK}
\alias{pUmbrPK}
\title{
Function to compute the P-value for the observed Mack-Wolfe Peak Known A_p distribution.
}
\description{
The function generalizes Harding's (1984) algorithm to quickly generate the distribution of A_p.
}
\usage{
pUmbrPK(x,peak=NA,g=NA,method=NA, n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing the data.}
 \item{peak}{An integer representing the known peak among the k data groups.}
 \item{g}{If x is a vector, g is a required vector of group labels. Otherwise, not used.}
 \item{method}{Either "Exact", "Monte Carlo", or "Asymptotic", indicating the desired distribution. When method=NA, and there are ties in the data, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used. When method=NA and there are no ties in the data, if sum(n)<=200, the "Exact" method will be used to compute the A_p distribution. Otherwise, the "Asymptotic" method will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the groups of data can be entered in either of two ways. For data a=1,2 and b=3,4,5 the following are equivalent:

\code{pUmbrPK(x=list(c(1,2),c(3,4,5)))}
\code{pUmbrPK(x=c(1,2,3,4,5),g=c(1,1,2,2,2))}

}
\value{
Returns a list with "NSM3Ch6p" class containing the following components:
\item{n}{a vector containing the number of observations in each of the data groups}
\item{obs.stat}{the observed A_p statistic}
\item{p.val}{the upper tail P-value}
}
\references{
Harding, E. F. "An efficient, minimal-storage procedure for calculating the Mann-Whitney U, generalized U and similar distributions." Applied statistics (1984): 1-6.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 6.3 Fasting Metabolic Rate of White-Tailed Deer
x<-c(36,33.6,26.9,35.8,30.1,31.2,35.3,39.9,29.1,43.4,44.6,54.4,48.2,55.7,50,53.8,53.9,62.5,46.6,
44.3,34.1,35.7,35.6,31.7,22.1,30.7)
g<-c(rep(1,7),rep(2,3),rep(3,5),rep(4,4),rep(5,4),rep(6,3))

pUmbrPK(x,4,g,"Exact")
pUmbrPK(x,4,g,"Asymptotic")
}
\keyword{Mack-Wolfe}
\keyword{Peak Known}