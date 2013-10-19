\name{multComb}
\alias{multComb}
\title{
Combinations of the first n integers in k groups
}
\description{
This is a function, used for generating the permutations used for the Exact distribution of many of the statistical procedures in Hollander, Wolfe, Chicken - Nonparametric Statistical Methods Third Edition, to generate possible combinations of the first n=n1+n2+...+nk integers within k groups.
}
\usage{
multComb(n.vec)
}
\arguments{
  \item{n.vec}{Contains the group sizes n1,n2,...,nk}
}
\details{
The computations involved get very time consuming very quickly, so be careful not to use it for too many large groups.
}
\value{
Returns a matrix of n!/(n1!*n2!*...*nk!) rows, where each row represents one possible combination.
}
\author{
Grant Schneider
}
\examples{
##What are the ways that we can group 1,2,3,4,5 into groups of 2, 2, and 1?
multComb(c(2,2,1))

##Another example, with four groups
multComb(c(2,2,3,2))
}
\keyword{Combinations}
\keyword{k groups}
