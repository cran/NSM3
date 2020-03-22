\name{MillerJack}
\alias{MillerJack}
\title{
Miller Jackknife
}
\description{
Function to compute the Miller Jackknife Q statistic.
}
\usage{
MillerJack(x,y=NULL)
}
\arguments{
  \item{x}{Either a vector containing the first group of data (X) or a matrix containing both groups of data.
}
\item{y}{If x is a vector, y is a vector containing the second group of data (Y). Otherwise, not used.}
}
\value{
Returns the observed Q statistic.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 5.2 Southern Armyworm and Pokeweed
kentucky.pokeweed<-c(6.2,5.9,8.9,6.5,8.6)
florida.pokeweed<-c(9.5,9.8,9.5,9.6,10.3)
MillerJack(kentucky.pokeweed,florida.pokeweed)
}
\keyword{Miller}
\keyword{Jackknife}