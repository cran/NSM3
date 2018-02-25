\name{theil}
\alias{theil}
\title{
Function to estimate and perform tests on the slope and intercept of a simple linear model.
}
\description{
This code estimates and performs tests on the slope and intercept of a simple linear model.  Based on chapter 9 of Hollander, Wolfe & Chicken, Nonparametric Statistical Methods, 3e.
}
\usage{
theil(x=NULL, y=NULL, alpha=0.05, beta.0=0, type="t", 
      example=FALSE, r=3, slopes=F, doplot=TRUE) 
}
\arguments{
 \item{x}{first data vector}
 \item{y}{second data vector}
 \item{alpha}{the significance level}
 \item{beta.0}{the null hypothesized value}
 \item{type}{can be "t" (two-sided), "u" (upper) or "l" (lower).  The type refers both to the test and the confidence interval.}
 \item{example}{if true, will analyze the data from Example 9.1}
 \item{r}{the number of places for rounding.  Increase it if your P-values are coming out as 0 or 1.}
 \item{slopes}{if true, will print all n(n-1)/2 slopes}
 \item{doplot}{if true, will plot the data and estimated line}
}
\value{
Returns a list with "NSM3Ch9ChickFn" class containing the following components:

\item{alpha}{same as input argument}
\item{beta.0}{same as input argument}
\item{type}{same as input argument}
\item{r}{same as input argument}
\item{slopes}{same as input argument}
\item{C.stat}{the observed C statistic}
\item{C.bar}{the observed C.bar statistic}
\item{alpha.hat}{the observed alpha.hat statistic}
\item{beta.hat}{the observed beta.hat statistic}
\item{slopes.table}{table containing all n(n-1)/2}
\item{p.val}{the P-value corresponding to the selected type of test/confidence interval}
\item{L}{the lower endpoint of the confidence interval}
\item{U}{the upper endpoint of the confidence interval}

}

\author{
Eric Chicken
}
\examples{
##Example 9.1 Hollander-Wolfe-Chicken##
theil (x, y, example=TRUE, slopes=TRUE)
}
\keyword{theil}
