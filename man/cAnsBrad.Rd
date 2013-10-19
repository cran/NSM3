\name{cAnsBrad}
\alias{cAnsBrad}
\title{
Function to compute a critical value for the Ansari-Bradley C distribution.
}
\description{
This function uses pAnsari and qAnsari from the base stats package to compute the critical value for the Ansari-Bradley C distribution at (or typically in the "Exact" case, close to) the given alpha level.  The program is reasonably quick for large data, well after the asymptotic approximation suffices, so Monte Carlo methods are not included.
}
\usage{
cAnsBrad(alpha, m, n, method = NA, n.mc = 10000)
}
\arguments{
  \item{alpha}{
  A numeric value between 0 and 1.
}
  \item{m}{
  A numeric value indicating the size of the first data group (X).
}
  \item{n}{
  A numeric value indicating the size of the second data group (Y).
}
  \item{method}{
  Either "Exact" or "Asymptotic", indicating the desired distribution. When method=NA, if m+n<=200, the "Exact" method will be used to compute the C distribution. Otherwise, the "Asymptotic" method will be used.
}
  \item{n.mc}{
  Not used.  Only included for standardization with other critical value procedures in the NSM3 package.
}
}
\value{
Returns a list with "NSM3Ch5c" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y)}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U (if method="Exact")}
\item{cutoff.L}{lower tail cutoff at or below user-specified alpha}
\item{true.alpha.L}{true alpha level corresponding to cutoff.L (if method="Exact")}
}
\author{
Grant Schneider
}
\seealso{
Also see \code{ansari.test()}
}

\references{
This function uses the source code ansari.c from the stats package by: R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.

}

\examples{
##Hollander, Wolfe, Chicken - NSM3 - Example 5.1 (Serum Iron Determination):
cAnsBrad(0.05,20,20,"Asymptotic")
cAnsBrad(0.05,20,20,"Exact")

##Bigger data
cAnsBrad(0.05,100,100,"Exact")
}
\keyword{Ansari-Bradley}