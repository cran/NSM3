\name{pFligPoli}
\alias{pFligPoli}
\title{
Fligner-Policello
}
\description{
Function to compute the P-value for the observed Fligner-Policello U statistic.
}
\usage{
pFligPoli(x,y=NA,g=NA,method=NA,n.mc=10000)
}
\arguments{
\item{x}{Either a list or a vector containing either all or the first group of data.}
\item{y}{If x contains the first group of data, y contains the second group of data. Otherwise, not used.}
\item{g}{If x contains a vector of all of the data, g is a vector of 1's and 2's corresponding to group labels. Otherwise, not used.}
\item{method}{Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the two groups of data can be entered in any of three ways. For data a=1,2 and b=3,4 all of the following are equivalent:

\code{pFligPoli(x=c(1,2),y=c(3,4))}
\code{pFligPoli(x=list(c(1,2),c(3,4)))}
\code{pFligPoli(x=c(1,2,3,4),g=c(1,1,2,2))}

}
\value{
Returns a list with "NSM3Ch5p" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{obs.stat}{the observed U statistic}
\item{p.val}{upper tail P-value}
\item{two.sided}{two-sided P-value}
}
\author{
Grant Schneider
}
\examples{
##Hollander, Wolfe, Chicken Example 4.5 Plasma Glucose in Geese
plasma.glucose<-list(healthy.geese = c(297, 340, 325, 227, 277, 337, 
250, 290), poisoned.geese = c(293, 291, 289, 430, 510, 353, 318
))

pFligPoli(plasma.glucose)

}
\keyword{Fligner-Policello}