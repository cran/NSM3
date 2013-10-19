\name{pHoeff}
\alias{pHoeff}
\title{Hoeffding's D}
\description{
Function to approximate the distribution of Hoeffding's D statistic using a Monte Carlo Sample under the null hypothesis.  This code follows section 8.6 of Hollander, Wolfe & Chicken, Nonparametric Statistical Methods, 3e. This calls HoeffD, a small bit of code that produces the value of D without any inference.  It is intended for small sample sizes n only.  For large n, use the asymptotic equivalence of D to the Blum-Kliefer-Rosenblatt statistic in the R package "Hmisc", command "hoeffd".
}
\usage{
pHoeff(n=5, reps=10000, r=4)
}
\arguments{
  \item{n}{the sample size}
  \item{reps}{the number of Monte Carlo runs to produce}
  \item{r}{the number of digits for rounding the results}
}
\value{
Returns a matrix containing the Monte Carlo distribution of the D statistic.
}
\author{
Eric Chicken
}
\seealso{
Also see the Hmisc package.
}
\examples{
pHoeff(n=5, reps=10000, r=4)
pHoeff(n=10, reps=1000, r=5)
}
\keyword{Hoeffding}
\keyword{Monte Carlo}