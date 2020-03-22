\name{HollBivSym}
\alias{HollBivSym}
\title{
Hollander Bivariate Symmetry
}
\description{
Function to compute the Hollander A statistic for testing bivariate symmetry.
}
\usage{
HollBivSym(x,y=NULL)
}
\arguments{
  \item{x}{Either a matrix containing both groups of data or a vector containing the first group of data.
}
\item{y}{If x is a vector, y is a required vector containing the second group of data. Otherwise, not used.
}
}
\details{
The data entry is intended to be flexible, so that the data can be entered in either of two ways. The following are equivalent:

\code{HollBivSym(x=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T))}
\code{HollBivSym(x=c(1,3,5),y=c(2,4,6))}
}
\value{
Returns the observed Hollander A statistic.
}
\author{
Grant Schneider
}
\examples{

##Hollander-Wolfe-Chicken Table 3.16 example
recipient<-c(61.4,63.3,63.7,80,77.3,84,105)
donor<-c(70.8,89.2,65.8,67.1,87.3,85.1,88.1)

HollBivSym(recipient,donor)

##Or, equivalently
table3.16<-matrix(c(61.4,63.3,63.7,80,77.3,84,105,70.8,89.2,65.8,67.1,87.3,85.1,88.1),ncol=2)
HollBivSym(table3.16)
}
\keyword{Hollander}