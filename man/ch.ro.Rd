\name{ch.ro}
\alias{ch.ro}
\title{
Campbell-Hollander
}
\description{
Function to compute the Campbell-Hollander estimator G-hat
}
\usage{
ch.ro (x,n,alpha,mu,...)
}
\arguments{
    \item{x}{a vector of data of length r}
    \item{n}{the sample size}
    \item{alpha}{the degrees of confidence in mu}
    \item{mu}{the prior guess of the unknown P (a pdf)}
    \item{...}{all of the arguments needed for mu}

}
\value{
\item{G.hat}{ estimate of the rank order G} 
}
\references{
See Section 16.3 of Hollander, Wolfe, Chicken - Nonparametric Statistical Methods 3.
}
\author{
Rachel Becvarik
}
\examples{
##Hollander-Wolfe-Chicken Example 16.2 Swimming in the Women's 50 yard Freestyle
freestyle<-c(22.43, 21.88, 22.39, 22.78, 22.65, 22.60)
ch.ro(freestyle,12,10,pnorm,22.52,.24)
}
\keyword{Campbell-Hollander}
