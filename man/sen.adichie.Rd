\name{sen.adichie}
\alias{sen.adichie}
\title{
Function to test for parallel lines.
}
\description{
This code tests for parallel lines based on chapter 9 of Hollander, Wolfe, & Chicken, Nonparametric Statistical Methods, 3e.
}
\usage{
sen.adichie(z, example=F, r=3) 
}
\arguments{
\item{z}{a list of paired vectors.  Each item in the list is a set of two paired vectors in the form of a matrix.  The first column of each matrix is the x vector, the second in the y vector.}
\item{example}{if true, analyzes the data from Example 9.5}
\item{r}{determines the amount of rounding. Increase it if your P-values are coming out as 0 or 1.}
}
\author{
Eric Chicken
}
\examples{
##Example 9.5 Hollander-Wolfe-Chicken##
sen.adichie(example=TRUE)
}
\keyword{Sen}
\keyword{Adichie}