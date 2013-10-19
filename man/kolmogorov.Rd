\name{kolmogorov}
\alias{kolmogorov}
\title{
Kolmogorov
}
\description{
Function to compute the asymptotic P-value for the observed Kolmogorov D statistic.
}
\usage{
kolmogorov(x,fnc,...)
}
\arguments{
\item{x}{a vector of data of length n}
\item{fnc}{the functional form of the pdf of F0.  The first argument must be the data.}
\item{...}{all the parameters besides the data that fnc needs to operate.  (See below for an example using pnorm and pexp)}
}
\value{
The function returns a list with two elements:
\item{D}{the value of the Kolmogorov statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
velocity<-c(12.8, 12.9, 13.3, 13.4, 13.7, 13.8, 14.5)
kolmogorov(velocity,pnorm, mean=14,sd=2)
kolmogorov(velocity,pexp,1/2)
}
\keyword{Kolmogorov}
\keyword{Asymptotic}