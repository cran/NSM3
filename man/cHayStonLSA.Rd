\name{cHayStonLSA}
\alias{cHayStonLSA}
\title{
Computes a critical value for the Hayter-Stone W* asymptotic distribution.
}
\description{
This function computes the critical value for the Hayter-Stone W* asymptotic distriburion at the given alpha level.
}
\usage{
cHayStonLSA(alpha,k,delta=.001)
}
\arguments{
   \item{alpha}{
  A numeric value between 0 and 1.
  }
  \item{k}{
  A numeric value indicating the number of the data groups (with assumed equal sizes).
}
  \item{delta}{
  Increment used to create the grid on which the distribution will be approximated.
  }

}
\details{
The Asymptotic distribution requires that all (unspecified) group sizes are equal.
}
\value{
Returns the cutoff (based on the specified grid) with upper tail probability nearest to alpha.
}
\references{
Hayter, Anthony J., and Wei Liu. "Exact calculations for the one-sided studentized range test for testing against a simple ordered alternative." Computational statistics & data analysis 22.1 (1996): 17-25.
}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 6.7 Motivational Effect of Knowledge of Performance:
cHayStonLSA(.0553,3,delta=0.01)

##Section preceding Example 6.7 (explaining LSA)
cHayStonLSA(.05,6,delta=0.01)
}
\keyword{Hayter-Stone LSA}