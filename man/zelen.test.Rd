\name{zelen.test}
\alias{zelen.test}
\title{
Function to perform Zelen's test.
}
\description{
Zelen's test based on section 10.4 of Hollander, Wolfe, & Chicken, Nonparametric Statistical Methods, 3e.
}
\usage{
zelen.test(z, example=F, r=3)
}
\arguments{
\item{z}{ data as an array of k 2x2 matrices.  Small data sets only!}
\item{example}{if true, analyzes the data from comment 24 of Chapter 10}
\item{r}{determines the amount of rounding. Increase it if your P-values are coming out as 0 or 1.}
}
\author{
Eric Chicken
}
\examples{
##Chapter 10 Coment 24 Hollander-Wolfe-Chicken##
zelen.test(example=TRUE)
}
\keyword{zelen}