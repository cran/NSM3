\name{RSS}
\alias{RSS}
\title{
Ranked-Set Sample
}
\description{
Function to obtain a ranked-set sample of given set size and number of cycles based on a specified auxiliary variable.
}
\usage{
RSS(k,m,ranker)
}
\arguments{
\item{k}{set size}
\item{m}{number of cycles}
\item{ranker}{auxiliary variable used for judgment ranking}
}
\value{
Returns a vector of the indices corresponding to the observations selected to be in the RSS.
}
\author{
Grant Schneider
}
\examples{
##Simulate 100 observations of a response variable we are interested in 
##and an auxiliary variable we use for ranking

set.seed(1)
response<-rnorm(100)
auxiliary<-rnorm(100)

##Get the indices for a ranked-set sample with set size 3 and 2 cycles
RSS(2,3,auxiliary) #Tells us to measure observations 2, 19, 32,..., 91

##Alternatively, get the responses for those observations. 
##In practice, response will not be available ahead of time.
response[RSS(2,3,auxiliary)]
}
\keyword{Ranked-Set Sample}
\keyword{RSS}