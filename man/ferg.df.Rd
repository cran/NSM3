\name{ferg.df}
\alias{ferg.df}
\title{
Ferguson's Estimator
}
\description{
Function to compute an approximation of Ferguson's estimator mu_n.
}
\usage{
ferg.df(x, alpha, mu, npoints,...)
}
\arguments{
  \item{x}{a vector of data of length n}
  \item{alpha}{the degree of confidence in mu}
  \item{mu}{the prior guess of the unknown P (a pdf)}
  \item{npoints}{the number of estimated points returned}
  \item{...}{all of the arguments needed for mu}

}
\value{
The function returns a vector of length num.points for Ferguson's estimator. 
}
\references{
See Section 16.2 of Hollander, Wolfe, Chicken - Nonparametric Statistical Methods 3.
}
\author{
Rachel Becvarik
}
\examples{
##Hollander-Wolfe-Chicken Figure 16.2
framingham<-c(2273, 2710, 141, 4725, 5010, 6224, 4991, 458, 1587, 1435, 2565, 1863)
plot.ecdf(framingham)
lines(sort(framingham),pexp(sort(framingham), 1/2922), lty=3)
temp.x = seq(min(framingham), max(framingham), length.out=100)
lines(temp.x,ferg.df(sort(framingham), 4, npoints=100,pexp,1/2922), col=2, type="s", lty=2)
legend("bottomright",  lty=c(1,3,2), legend=c("ecdf", "prior", "ferguson"), col=c(1,1,2))
}
\keyword{Ferguson}