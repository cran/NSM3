\name{pKW}
\alias{pKW}
\title{
Kruskal-Wallis
}
\description{
Function to compute the P-value for the observed Kruskal-Wallis H statistic.
}
\usage{
pKW(x,g=NA, method=NA, n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing the data.}
 \item{g}{If x is a vector, g is a required vector of group labels. Otherwise, not used.}
  \item{method}{Either "Exact", "Monte Carlo", or "Asymptotic", indicating the desired distribution. When method=NA and ties are not present, "Exact" will be used. When method=NA and ties are present, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the groups of data can be entered in either of two ways. For data a=1,2 and b=3,4,5 the following are equivalent:

\code{pKW(x=list(c(1,2),c(3,4,5)))}
\code{pKW(x=c(1,2,3,4,5),g=c(1,1,2,2,2))}

}
\value{
Returns a list with "NSM3Ch6p" class containing the following components:
\item{n}{a vector containing the number of observations in each of the data groups}
\item{obs.stat}{the observed H statistic}
\item{p.val}{upper tail P-value}
}
\author{
Grant Schneider
}
\seealso{
Also see \code{kruskal.test()}.
}
\examples{
##Hollander-Wolfe-Chicken Example 6.1 Half-Time of Mucociliary Clearance
mucociliary<-list(Normal = c(2.9, 3, 2.5, 2.6, 3.2), Obstructive = c(3.8, 
2.7, 4, 2.4), Asbestosis = c(2.8, 3.4, 3.7, 2.2, 2))

pKW(mucociliary)

}
\keyword{Kruskal-Wallis}
