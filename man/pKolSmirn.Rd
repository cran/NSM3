\name{pKolSmirn}
\alias{pKolSmirn}
\title{
Function to copute the P-value for the observed Kolmogorov-Smirnov J statistic.
}
\description{
This function uses psmirnov2x from the base stats package to compute the J statistic and corresponding P-value. The program is reasonably quick for large data, well after the asymptotic approximation suffices, so Monte Carlo methods are not included. This function primarily serves as a wrapper to the ks.test function with the output standardized to the format of the other functions included in the NSM3 package.
}
\usage{
pKolSmirn(x,y=NA,g=NA,method=NA,n.mc=10000)
}
\arguments{
 \item{x}{Either a list or a vector containing either all or the first group of data.}
\item{y}{If x contains the first group of data, y contains the second group of data. Otherwise, not used.}
\item{g}{If x contains a vector of all of the data, g is a vector of 1's and 2's corresponding to group labels. Otherwise, not used.}
\item{method}{Either "Exact" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used.  
}
\item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the two groups of data can be entered in any of three ways. For data a=1,2 and b=3,4 all of the following are equivalent:

\code{pKolSmirn(x=c(1,2),y=c(3,4))}
\code{pKolSmirn(x=list(c(1,2),c(3,4)))}
\code{pKolSmirn(x=c(1,2,3,4),g=c(1,1,2,2))}

}
\value{
Returns a list with "NSM3Ch5p" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{obs.stat}{the observed C statistic}
\item{p.val}{upper tail P-value}
}
\author{
Grant Schneider
}
\seealso{
Also see \code{ks.test()}.
}
\examples{
##Hollander-Wolfe-Chicken Example 5.4 Effect of Feedback on Salivation Rate:
feedback<-c(-0.15, 8.6, 5, 3.71, 4.29, 7.74, 2.48, 3.25, -1.15, 8.38)
no.feedback<-c(2.55, 12.07, 0.46, 0.35, 2.69, -0.94, 1.73, 0.73, -0.35, -0.37)
pKolSmirn(x=feedback,y=no.feedback)
}
\keyword{Kolmogov-Smirnov}

