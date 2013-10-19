\name{mrl}
\alias{mrl}
\title{
Mean Residual Life
}
\description{
Function to return the mean residual life along with Hall and Wellner's upper and lower bounds.
}
\usage{
mrl(data, alpha, main=NULL, ylim=NULL, xlab=NULL,...)
}
\arguments{
  \item{data}{a vector of survival times}
\item{alpha}{(1-alpha) is the approximate coverage probability for the confidence band.}
\item{main}{title of the plot.  The default is "Plot of Mean Residual Life and bounds".}
\item{ylim}{the limits of the y-axis.  The default is to include all points in the plotting range.}
\item{xlab}{the label for the x-axis.  The default is Time.}
\item{...}{additional plotting options}
}
\value{
The function returns a list with three vectors:
\item{PM}{the mean residual life}
\item{PMU}{upper bound for the mean residual life}
\item{PML}{lower bound for the mean residual life}
}
\author{
Rachel Becvarik
}
\examples{
leukemia<-c(7, 429, 579, 968, 1877, 47, 440, 581, 1077, 1886, 58,
445,  650, 1109, 2045, 74, 455, 702, 1314, 2056, 177, 468,
715, 1334, 2260, 232, 495, 779, 1367, 2429, 273, 497, 881,
1534, 2509, 285, 532, 900, 1712, 317,  571, 930, 1784)
mrl(leukemia, .05)
}
\keyword{Hall and Wellner}
\keyword{MRL}