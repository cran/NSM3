\name{newbet}
\alias{newbet}
\title{
Hollander-Proschan T*
}
\description{
Function to compute the asymptotic P-value for the observed Hollander-Proschan T* statistic.
}
\usage{
newbet(x)
}
\arguments{
  \item{x}{a vector of data of length n}
}
\value{
The function returns a list with two elements:
\item{T}{the value of the Hollander-Proschan statistic}
\item{T.star}{the standardized value of the Hollander-Proschan statistic}
\item{p}{the corresponding probability}
}
\author{
Rachel Becvarik
}
\examples{
table11.4<-c(194,15,41,29,33,181)
newbet(table11.4)
}
\keyword{Hollander-Proschan}

