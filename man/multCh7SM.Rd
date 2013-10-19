\name{multCh7SM}
\alias{multCh7SM}
\title{
Possible arrangements by row a matrix, where NA values are ignored
}
\description{
Similar to \code{multCh7}, this function will generate all of the possible arrangements of the data by row within a matrix, except for NA values, which will remain fixed. This function is used in pSkilMack and cSkilMack to generate the Exact distribution. For a given matrix of with k1,...kn non-missing values, this will give k1!*k2!*...*kn! possible arrangements
}
\usage{
multCh7SM(our.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{our.matrix}{The matrix containing the data (including NA values) which will be rearranged by row.}
}
\details{
The computations involved get very time consuming very quickly, so be careful not to use it for too large of a matrix.
}
\value{Returns an array, containing k1!*k2!*...*kn! distinct matrices of the same size as \code{our.matrix}}
\author{
Grant Schneider
}
\examples{
##Get a matrix with some NA's
our.matrix<-matrix(c(NA,1,2,3,5,7,NA,NA,11),ncol=3,byrow=TRUE)
##Get every possible arrangement by row, treating the NA's as fixed
multCh7SM(our.matrix)
}
\keyword{Row arrangements NA fixed}