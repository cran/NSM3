\name{pBohnWolfe}
\alias{pBohnWolfe}
\title{
Function to compute the P-value for the observed Bohn-Wolfe U statistic.
}
\description{
This function computes the U statistic and then uses Monte Carlo sampling to compute the corresponding P-value. The Monte Carlo samples are simulated based on the order statistics of a uniform(0,1) distribution.
}
\usage{
pBohnWolfe(x,y,k,q,c,d,method="Monte Carlo",n.mc=10000)
}
\arguments{
  \item{x}{A vector containing the data in the first group.}
  \item{y}{A vector containing the data in the Second group.}
  \item{k}{
  A numeric value indicating the set size of the first data group in the RSS (X).
}
  \item{q}{
  A numeric value indicating the set size of the second data group in the RSS (Y).
}
  \item{c}{
  A numeric value indicating the number of cycles for the first data group in the RSS (X).
}
  \item{d}{
  A numeric value indicating the number of cycles for the second data group in the RSS (Y).
}
  \item{method}{
For this procedure, method is currently set automatically to "Monte Carlo" as the only option that is available.  For standardization with other critical value procedures in the NSM3 package, "Asymptotic" and "Exact" will be supported in future versions.
}
  \item{n.mc}{
  Number of Monte Carlo samples used to estimate the distribution of U.
}
}
\value{
Returns a list with "NSM3Ch5p" class containing the following components:
\item{m}{number of observations in RSS for the first data group (X)}
\item{n}{number of observations in RSS for the second data group (Y)}
\item{obs.stat}{the observed U statistic}
\item{p.val}{upper tail P-value}
}
\references{
Bohn, Lora L., and Douglas A. Wolfe. "Nonparametric two-sample procedures for ranked-set samples data." Journal of the American Statistical Association 87.418 (1992): 552-561
}
\author{
Grant Schneider
}
\examples{
##Hollander, Wolfe, Chicken Example 15.4 Body Mass Index:
male<-c(18.0, 20.5, 21.3, 21.3, 22.3, 23.8, 23.8, 24.6, 25.0, 25.2, 25.3, 25.9, 26.1, 27.0,
27.4, 27.4, 28.4, 29.4, 29.6, 32.8)
female<-c(17.2, 17.8, 19.9, 20.0, 21.7, 22.0, 22.3, 23.1, 23.9, 25.8, 27.1, 29.6, 30.1, 30.3,
30.7, 31.1, 35.2, 35.6, 38.1, 42.5)

pBohnWolfe(male,female,4,4,5,5)
##To use more Monte Carlo samples:
#pBohnWolfe(male,female,4,4,5,5,n.mc=100000)

}
\keyword{Bohn-Wolfe}
\keyword{Ranked-Set Sample}