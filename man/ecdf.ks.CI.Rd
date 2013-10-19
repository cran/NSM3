\name{ecdf.ks.CI}
\alias{ecdf.ks.CI}
\title{
Kolmogorov's Confidence Band
}
\description{
Function to compute and plot Kolmogorov's 95\% confidence band for the distribution function F(x). This code is adapted from the code by Kjetil Halvorsen found at: https://stat.ethz.ch/pipermail/r-help/2003-July/036643.html
}
\usage{
ecdf.ks.CI(x, main = NULL, sub = NULL, xlab = deparse(substitute(x)), ...)
}
\arguments{
\item{x}{a vector of data of length n}
\item{main}{the title of the plot.  The default is ecdf(x) + 95\% K.S.Bands}
\item{sub}{subtitle, as used in the function plot()}
\item{xlab}{the label for the x-axis of the plot.  The default is x.}
\item{...}{any additional plotting options}
}
\value{
The function returns a list with three elements:
\item{lower}{the values of the lower part of the confidence band}
\item{upper}{the values of the upper part of the confidence band}
\item{D}{the value of Kolmogorov's D statistic}  
}
\author{
Rachel Becvarik
}
\note{
This function also plots the confidence bands.
}
\examples{
methyl<-c(42, 43, 51, 61, 66, 69, 71, 81, 82, 82)
ecdf.ks.CI(methyl)

ecdf.ks.CI(methyl, lwd=2, main="KS Confidence Bands") 

}
\keyword{Kolmogorov}
\keyword{CI}