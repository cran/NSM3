\name{NSM3-package}
\alias{NSM3-package}
\alias{NSM3}
\docType{package}
\title{
This package implements nonparametric functions to accompany Hollander, Wolfe, and Chicken - Nonparametric Statistical Methods 3rd edition.
}
\description{
Generally, the function are of two forms, "c" and "p". The "c" functions will compute a cut-off for a given alpha level, while the "p" functions will compute the test statistic and P-value for a given data set. Typically, Exact, Monte Carlo, and Asymptotic methods are available, but this may vary on a case by case basis. The functions are able to handle ties, a significant extension to the methods available previously.
}
\details{
\tabular{ll}{
Package: \tab NSM3\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-02\cr
License: \tab GPL 2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Grant Schneider, Eric Chicken, and Rachel Becvarik

Maintainer: Grant Schneider <schneider.393@osu.edu>

}
\references{
Hollander, Wolfe, and Chicken - Nonparametric Statistical Methods, Third Edition
}
\keyword{Nonparametric}
\keyword{Hollander}
\keyword{Wolfe}
\keyword{Chicken}
\examples{
##Hollander, Wolfe, Chicken - NSM3 - Example 5.1 (Serum Iron Determination):
cAnsBrad(0.05,20,20,"Asymptotic")
cAnsBrad(0.05,20,20,"Exact")

##Bigger data
cAnsBrad(0.05,100,100,"Exact")

}
