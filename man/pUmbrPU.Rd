\name{pUmbrPU}
\alias{pUmbrPU}
\title{
Mack-Wolfe Peak Unknown
}
\description{
Function to compute the P-value for the observed Mack-Wolfe Peak Unknown A_p-hat distribution.
}
\usage{
pUmbrPU(x,g=NA,method=NA, n.mc=10000)
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

\code{pUmbrPU(x=list(c(1,2),c(3,4,5)))}
\code{pUmbrPU(x=c(1,2,3,4,5),g=c(1,1,2,2,2))}

}
\value{
Returns a list with "NSM3Ch6p" class containing the following components:
\item{n}{a vector containing the number of observations in each of the data groups}
\item{obs.stat}{the observed A_p-hat statistic}
\item{p.val}{the upper tail P-value}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 6.4 Learning Comprehension and Age
wechsler<-list("16-19"=c(8.62,9.94,10.06),"20-34"=c(9.85,10.43,11.31),"35-54"=c(9.98,10.69,11.40),
"55-69"=c(9.12,9.89,10.57),"70+"=c(4.80,9.18,9.27))

#pUmbrPU(wechsler,method="Monte Carlo",n.mc=20000)
pUmbrPU(wechsler,method="Monte Carlo",n.mc=1000)
}
\keyword{Mack-Wolfe}
\keyword{Peak Unknown}