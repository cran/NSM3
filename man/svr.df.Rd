\name{svr.df}
\alias{svr.df}
\title{
Susarla-van Ryzin
}
\description{
Function to compute the Susarla-van Ryzin estimator 
}
\usage{
svr.df (z, delta, lambda.hat=0.001, alpha = 3, npoints=2053)
}
\arguments{
\item{z}{the vector of zi = min(Xi, Yi)}
\item{delta}{the vector of indicators which is 1 when Xi<=Yi and 0 otherwise}
\item{lambda.hat}{the estimate of lambda from the data}
\item{alpha}{the degree of faith in F0}
\item{npoints}{the number of estimated points returned}

}
\value{
Returns a list containing:
\item{x}{the x values}
\item{F.hat}{the Susarla-van Ryzin estimator}
}
\author{
Rachel Becvarik
}
\note{
Requires the survival library.
}
\examples{
hodgkins.affected<-matrix(c(1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1,0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
0, 1, 346, 141, 296, 1953, 1375, 822, 2052, 836, 1910, 419,  107, 570, 312,1818, 364, 401, 1645,
330, 1540, 688, 1309, 505, 1378, 1446, 86),nrow=2,byrow=TRUE)
svr.df(hodgkins.affected[2,], hodgkins.affected[1,])
}
\keyword{Susarla}
\keyword{van Ryzin}
