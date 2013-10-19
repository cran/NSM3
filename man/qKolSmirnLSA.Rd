\name{qKolSmirnLSA}
\alias{qKolSmirnLSA}
\title{
Quantile function for the asymptotic distribution of the Kolmogorov-Smirnov J* statistic.
}
\description{
This function computes the Q() function defined in Section 5.4 of Hollander, Wolfe, and Chicken on a grid and then searches for the cutoff based on alpha.
}
\usage{
qKolSmirnLSA(alpha)
}
\arguments{
  \item{alpha}{
      A numeric value between 0 and 1.
  }
}
\value{
Returns the upper tail cutoff at or below user-specified alpha
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Section 5.4 LSA
qKolSmirnLSA(.05)
}
\keyword{Kolmogorv-Smirnov}
\keyword{Asymptotic}