\name{owa}
\alias{owa}
\title{
Ordered Walsh Averages
}
\description{
Function to compute the ordered Walsh averages and the value of the Hodges-Lehmann estimator

}
\usage{
owa(x,y)
}
\arguments{
  \item{x}{ first vector of data of length n}
  \item{y}{ second vector of data of length n}
}
\value{
Returns a list containing:
\item{owa}{the ordered Walsh averages}
\item{h.l}{the value of the Hodges-Lehmann estimator}
}
\author{
Rachel Becvarik
}
\examples{
##Hollander-Wolfe-Chicken Example 3.3 
x<-c(1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y<-c(0.878, 0.647, 0.598, 2.050, 1.060, 1.290, 1.060, 3.140, 1.290)
owa(x,y)
}
\keyword{Ordered Walsh Averages}
\keyword{Hodges-Lehmann}