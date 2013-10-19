\name{pPairedWilcoxon}
\alias{pPairedWilcoxon}
\title{
Paired Wilcoxon
}
\description{
Function to extend wilcox.test to compute the (exact or Monte Carlo) P-value for paired Wilcoxon data in the presence of ties.
}
\usage{
pPairedWilcoxon(x,y=NA,g=NA,method=NA,n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing either all or the first group of data.}
\item{y}{If x contains the first group of data, y contains the second group of data. Otherwise, not used.}
\item{g}{If x contains a vector of all of the data, g is a vector of 1's and 2's corresponding to group labels. Otherwise, not used.}
  \item{method}{ Either "Exact" or "Monte Carlo", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the two groups of data can be entered in any of three ways. For data a=1,2 and b=3,4 all of the following are equivalent:

\code{pPairedWilcoxon(x=c(1,2),y=c(3,4))}
\code{pPairedWilcoxon(x=list(c(1,2),c(3,4)))}
\code{pPairedWilcoxon(x=c(1,2,3,4),g=c(1,1,2,2))}

}
\value{
Returns a list with "NSM3Ch5p" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{obs.stat}{the observed T+ statistic}
\item{p.val}{upper tail P-value}
}
\author{
Grant Schneider
}
\note{
If there are 0s in the Z values (the difference between X and Y), these will be removed and the calculations will be done based on the smaller sample size, as detailed section 3.1 of Hollander, Wolfe, and Chicken - NSM3.
}
\seealso{
Also see \code{stats::wilcox.test()}
}
\examples{
##Hollander-Wolfe-Chicken Example 3.1 Hamilton Depression Scale Factor IV
x <-c(1.83, .50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <-c(0.878, .647, .598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

wilcox.test(y,x,paired=TRUE,alternative="less")
pPairedWilcoxon(x,y)
}
\keyword{Wilcoxon}
\keyword{Paired}