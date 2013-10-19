\name{cHollBivSym}
\alias{cHollBivSym}
\title{
Hollander Bivariate Symmetry
}
\description{
Quantile function for the Hollander A distribution.
}
\usage{
cHollBivSym(alpha,d.mat,method=NA, n.mc=10000)
}
\arguments{
    \item{alpha}{
  A numeric value between 0 and 1.
}
\item{d.mat}{
The d matrix, explained below.}

\item{method}{ Either "Exact", "Monte Carlo" or "Asymptotic", indicating the desired distribution. When method=NA, "Exact" will be used if the number of permutations is 10,000 or less. Otherwise, "Monte Carlo" will be used. As Kepner and Randles (1984) and Hilton and Gee (1997) have found the large sample approximation to perform poorly, method="Asymptotic" will be treated as method=NA.
}
  \item{n.mc}{
  If method="Monte Carlo", the number of Monte Carlo samples used to estimate the distribution. Otherwise, not used.
}

}
\details{
The d matrix, d.mat, will be an n*n matrix of ones and zeroes, where the (i,j)th element is 1 if min(Xj,Yj)<max(Xi,Yi)<=max(Xj,Yj) and min(Xi,Yi)<=min(Xj,Yj), 0 otherwise. An illustration may be found in the example section of this document and Section 3.10 of Hollander, Wolfe, and Chicken - NSM3.
}
\value{
Returns a list with "NSM3Ch5c" class containing the following components:
\item{m}{number of observations in the first data group (X)}
\item{n}{number of observations in the second data group (Y) (equal to m, but included for standardization with other procedures)}
\item{cutoff.U}{upper tail cutoff at or below user-specified alpha}
\item{true.alpha.U}{true alpha level corresponding to cutoff.U}
}
\references{
Kepner, James L., and Ronald H. Randies. "Comparison of tests for bivariate symmetry versus location and/or scale alternatives." Communications in Statistics-Theory and Methods 13.8 (1984): 915-930.

Hilton, Joan F., and Lauren Gee. "The size and power of the exact bivariate symmetry test." Computational statistics & data analysis 26.1 (1997): 53-69.

}
\author{
Grant Schneider
}
\examples{
##Hollander-Wolfe-Chicken Example 3.11 Insulin Clearance in Kidney Transplants
x<-c(61.4,63.3,63.7,80,77.3,84,105)
y<-c(70.8,89.2,65.8,67.1,87.3,85.1,88.1)
obs.data<-cbind(x,y)
a.vec<-apply(obs.data,1,min)
b.vec<-apply(obs.data,1,max)
test<-function(r,c) {as.numeric((a.vec[c]<b.vec[r])&&(b.vec[r]<=b.vec[c])&&(a.vec[r]<=a.vec[c]))}
myVecFun <- Vectorize(test,vectorize.args = c('r','c')) 

d.mat<-outer(1:length(x), 1:length(x), FUN=myVecFun) 

##Cutoff based on the exact distribution
cHollBivSym(.10,d.mat)

}
\keyword{Hollander}
\keyword{Bivariate Symmetry}