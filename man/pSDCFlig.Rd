\name{pSDCFlig}
\alias{pSDCFlig}
\title{
Dwass, Steel, Critchlow, Fligner
}
\description{
Function to compute the P-value for the observed Dwass, Steel, Critchlow, Fligner W statistic.
}
\usage{
pSDCFlig(x,g=NA,method=NA,n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing the data.}
 \item{g}{If x is a vector, g is a required vector of group labels. Otherwise, not used.}
 \item{method}{Either "Exact", "Monte Carlo", or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the groups of data can be entered in either of two ways. For data a=1,2 and b=3,4,5 the following are equivalent:

\code{pSDCFlig(x=list(c(1,2),c(3,4,5)))}
\code{pSDCFlig(x=c(1,2,3,4,5),g=c(1,1,2,2,2))}

}
\value{
Returns a list with "NSM3Ch6MCp" class containing the following components:
\item{n}{a vector containing the number of observations in each of the k data groups}
\item{obs.stat}{the observed W statistic for each of the k*(k-1)/2 comparisons}
\item{p.val}{upper tail P-value corresponding to each W statistic}
}
\author{
Grant Schneider
}
\examples{
gizzards<-list(site.I=c(46,28,46,37,32,41,42,45,38,44),
              site.II=c(42,60,32,42,45,58,27,51,42,52),
              site.III=c(38,33,26,25,28,28,26,27,27,27),
              site.IV=c(31,30,27,29,30,25,25,24,27,30))
##Takes a little while 
#pSDCFlig(gizzards,method="Monte Carlo")

##Shorter version for demonstration
pSDCFlig(gizzards[1:2],method="Asymptotic")
}
\keyword{Dwass}
\keyword{Steel}
\keyword{Critchlow-Fligner}