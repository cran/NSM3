\name{RFPW}
\alias{RFPW}
\title{
Randles-Fligner-Policello-Wolfe
}
\description{
Function to compute the P-value for the observed Randles-Fligner-Policello-Wolfe V statistic.
}
\usage{
RFPW(z)
}
\arguments{
  \item{z}{
A vector containing the data.
}
}
\value{
Returns a list containing:
\item{obs.stat}{the observed V statistic}
\item{p.val}{the asymptotic two-sided P-value}
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 3.10 Percentage Chromium in Stainless Steel
table3.9.subset<-c(17.4,17.9,17.6,18.1,17.6)
RFPW(table3.9.subset)
}
\keyword{Randles-Fligner-Policello-Wolfe}
\keyword{Davis-Quade}