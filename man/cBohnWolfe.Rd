\name{cBohnWolfe}
\alias{cBohnWolfe}
\title{Function to compute a critical value for the Bohn-Wolfe U distribution.}
\description{
This function uses Monte Carlo sampling to compute the critical value for the Bohn-Wolfe U distribution at (or close to) the given alpha level. The Monte Carlo samples are simulated based on the order statistics of a uniform(0,1) distribution. 
}
\usage{
cBohnWolfe(alpha,k,q,c,d,method="Monte Carlo",n.mc=10000)
}
\arguments{
    \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{k}{
  A numeric value indicating the set size of the first data group in the RSS (X).
}
  \item{q}{
  A numeric value indicating the set size of the second data group in the RSS (Y).
}
  \item{c}{
  A numeric value indicating the number of cycles for the first data group in the RSS (X).
}
  \item{d}{
  A numeric value indicating the number of cycles for the second data group in the RSS (Y).
}
  \item{method}{
For this procedure, method is currently set automatically to "Monte Carlo" as the only option that is available.  For standardization with other critical value procedures in the NSM3 package, "Asymptotic" and "Exact" will be supported in future versions.
}
  \item{n.mc}{
  Number of Monte Carlo samples used to estimate the distribution of U.
}
}
\value{
Returns a list with "NSM3Ch5c" class containing the following components:
\item{m}{number of observations in RSS for the first data group (X)}
\item{n}{number of observations in RSS for the second data group (Y)}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U}
}
\references{
Bohn, Lora L., and Douglas A. Wolfe. "Nonparametric two-sample procedures for ranked-set samples data." Journal of the American Statistical Association 87.418 (1992): 552-561.
}
\author{
Grant Schneider
}
\examples{
cBohnWolfe(.0515,4,4,5,5)
cBohnWolfe(.0303,2,3,3,3)
}
\keyword{Bohn-Wolfe}
\keyword{RSS}