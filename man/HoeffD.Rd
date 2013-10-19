\name{HoeffD}
\alias{HoeffD}
\title{
Function to compute Hoeffding's D statistic for small sample sizes.
}
\description{
This will calculate Hoeffding's D statistic following section 8.6 of Hollander, Wolfe & Chicken, Nonparametric Statistical Methods, 3e.  Uses the correction for ties given at (8.92).
}
\usage{
HoeffD(x, y, example=FALSE) 
}
\arguments{
  \item{x}{first data vector}
  \item{y}{second data vector}
  \item{example}{if true, analyzes the data from Example 8.6}
}
\author{
Eric Chicken
}
\note{
This function is intended for small sample sizes n only.  For large n, use the asymptotic equivalence of D to the Blum-Kliefer-Rosenblatt statistic in the R package "Hmisc", command "hoeffd".
}
\examples{
##Example 8.6 Hollander-Wolfe-Chicken##
HoeffD(example=TRUE)
}
\keyword{Hoeffding}