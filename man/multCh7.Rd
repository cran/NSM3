\name{multCh7}
\alias{multCh7}
\title{
Possible arrangements by row for a matrix
}
\description{
Similar to \code{multComb}, this function will generate all of the possible arrangements of the data by row within a matrix. For a given matrix of n rows and k columns, this will give (k!)^n possible arrangements
}
\usage{
multCh7(our.matrix)
}
\arguments{
  \item{our.matrix}{The matrix containing the data which will be rearranged by row.}
}
\details{
The computations involved get very time consuming very quickly, so be careful not to use it for too large of a matrix.
}
\value{Returns an array, containing (k!)^n distinct matrices of the same size as \code{our.matrix}}
\author{
Grant Schneider
}
\note{
This function is used to generate the possible permutations for the Exact methods used in Chapter 7 of Hollander, Wolfe, and Chicken - Nonparametric Statistical Methods Third Edition.
}
\examples{
some.matrix<-matrix(c(1,2,7,4,5,9),ncol=3,byrow=TRUE)
multCh7(some.matrix)
}
\keyword{Row arrangements}