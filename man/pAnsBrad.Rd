\name{pAnsBrad}
\alias{pAnsBrad}
\title{
Function to compute the P-value for the observed Ansari-Bradley C statistic.
}
\description{
When there are no ties in the data, this function uses pansari and cansari from the base stats package to compute the C statistic and P-value ("Exact" or "Asymptotic").
The program is reasonably quick for large data in the absence of ties, well after the asymptotic approximation suffices, so Monte Carlo methods are not included.

When there are ties in the data, this function computes the C statistic and P-value ("Exact", "Monte Carlo", or "Asymptotic").

}
\usage{
pAnsBrad(x,y=NA,g=NA,method=NA,n.mc=10000)
}
\arguments{
\item{x}{Either a list or a vector containing either all or the first group of data.}
\item{y}{If x contains the first group of data, y contains the second group of data. Otherwise, not used.}
\item{g}{If x contains a vector of all of the data, g is a vector of 1's and 2's corresponding to group labels. Otherwise, not used.}
\item{method}{Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA and there are no ties in the data, "Exact" will be used.  When method=NA and there are ties in the data, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the two groups of data can be entered in any of three ways. For data a=1,2 and b=3,4 all of the following are equivalent:

\code{pAnsBrad(x=c(1,2),y=c(3,4))}
\code{pAnsBrad(x=list(c(1,2),c(3,4)))}
\code{pAnsBrad(x=c(1,2,3,4),g=c(1,1,2,2))}

}
\value{
Returns a list with "NSM3Ch5p" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{obs.stat}{the observed C statistic}
\item{p.val}{upper tail P-value}
\item{two.sided}{two-sided P-value}
}

\author{
Grant Schneider
}
\note{
If method="Monte Carlo" and there are no ties in the data, a warning is displayed and the "Exact" method is used.
}
\seealso{
Also see \code{ansari.test}.
}
\examples{
##Hollander, Wolfe, Chicken Example 5.1 Serum Iron Determination:
serum<-list(ramsay = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99, 101, 96, 97, 102, 107,
113, 116, 113, 110, 98),
jung.parekh = c(107, 108, 106, 98, 105, 103, 110, 105, 104, 100, 96, 108, 103, 104, 114, 114,
113, 108, 106, 99))


pAnsBrad(serum)

##or, equivalently:
pAnsBrad(serum$ramsay, serum$jung.parekh)

}
\keyword{Ansari-Bradley}