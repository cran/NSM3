\name{pHollBivSym}
\alias{pHollBivSym}
\title{
Hollander Bivariate Symmetry
}
\description{
Function to compute the P-value for the observed Hollander A statistic.
}
\usage{
pHollBivSym(x,y=NA,g=NA,method=NA,n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing either all or the first group of data.}
\item{y}{If x contains the first group of data, y contains the second group of data. Otherwise, not used.}
\item{g}{If x contains a vector of all of the data, g is a vector of 1's and 2's corresponding to group labels. Otherwise, not used.}
  \item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used. As Kepner and Randles (1984) and Hilton and Gee (1997) have found the large sample approximation to perform poorly, method="Asymptotic" will be treated as method=NA.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the two groups of data can be entered in any of three ways. For data a=1,2 and b=3,4 all of the following are equivalent:

\code{pHollBivSym(x=c(1,2),y=c(3,4))}
\code{pHollBivSym(x=list(c(1,2),c(3,4)))}
\code{pHollBivSym(x=c(1,2,3,4),g=c(1,1,2,2))}

}
\value{
Returns a list with "NSM3Ch5p" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{obs.stat}{the observed A statistic}
\item{p.val}{upper tail P-value}
}
\references{
Kepner, James L., and Ronald H. Randies. "Comparison of tests for bivariate symmetry versus location and/or scale alternatives." Communications in Statistics-Theory and Methods 13.8 (1984): 915-930.

Hilton, Joan F., and Lauren Gee. "The size and power of the exact bivariate symmetry test." Computational statistics & data analysis 26.1 (1997): 53-69.

}
\author{
Grant Schneider
}
\examples{

##Hollander-Wolfe-Chicken Example 3.11 Insulin Clearance in Kidney Transplants
x<-c(61.4,63.3,63.7,80,77.3,84,105)
y<-c(70.8,89.2,65.8,67.1,87.3,85.1,88.1)

##Exact p-value
pHollBivSym(x,y)

}
\keyword{Hollander}
\keyword{Bivariate Symmetry}