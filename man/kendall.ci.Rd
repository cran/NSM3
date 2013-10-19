\name{kendall.ci}
\alias{kendall.ci}
\title{
Function to produce a confidence interval for Kendall's tau.
}
\description{
Based on sections 8.3 and 8.4 of Hollander, Wolfe & Chicken, Nonparametric Statistical Methods, 3e.
}
\usage{
kendall.ci(x=NULL, y=NULL, alpha=0.05, type="t", bootstrap=F, B=1000, example=F) 
}
\arguments{
\item{x}{first data vector}
\item{y}{second data vector}
\item{alpha}{the significance level}
\item{type}{type of confidence interval. Can be "t" (two-sided), "u" (upper) or "l" (lower).}
\item{bootstrap}{if False, will find the asymptotic CI (as in section 8.3). If True, will find a bootstrap CI (as in section 8.4).}
\item{B}{the number of bootstrap replicates}
\item{example}{if True, will analyze data from Example 8.1}
}
\author{
Eric Chicken
}
\examples{
kendall.ci(example=TRUE)
}
\keyword{Kendall}
\keyword{CI}
\keyword{Confidence Interval}