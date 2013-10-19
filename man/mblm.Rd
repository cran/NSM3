\name{mblm}
\alias{mblm}
\title{
Fitting Median-Based Linear Models (from 'mblm' oackage)
}
\description{
This function is used to fit linear models based on Theil-Sen single median, or Siegel repeated medians.
}
\usage{
mblm(formula, dataframe, repeated = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{A formula of type y ~ x (only linear models are accepted)}
\item{dataframe}{Optional dataframe}
\item{repeated}{If set to true, model is computed using repeated medians. If false, a single median estimators are calculated}
}
\details{
This function is from the 'mblm' package, which is no longer available on CRAN.

Theil-Sen single median method computes slopes of lines crossing all possible pairs of points, when x coordinates differ. After calculating these n(n-1)/2 slopes (these value are true only if x is distinct), the median of them is taken as slope estimator. Next, the intercepts of n lines, crossing each point and having calculated slope are calculated. The median from them is intercept estimator.

Siegel repeated medians is more complicated. For each point, the slopes between it and the others are calcuated (resulting n-1 slopes) and the median is taken. This results in n medians and median from this medians is slope estimator. Intercept is calculated in similar way, for more information please take a look in function source.

The breakdown point of Theil-Sen method is about 29\%, Siegel extended it to 50\%, so these regression methods are very robust. Additionally, if the errors are normally distributed and no outliers are present, the estimators are very similar to classic least squares.
}
\value{
An object of class c("mblm","lm"), containing minimal set of data to perform basic operations, such as in case of lm model. Additionally, the return value contains 2 fields:

\item{slopes  }{The slopes (in single median), or medians of slopes (in repeated medians) between tested point pairs}
\item{intercepts}{The intercepts calculated}
}
\references{
Theil, H. (1950) A rank invariant method for linear and polynomial regression analysis. Nederl. Akad. Wetensch. Proc. Ser. A 53, 386-392 (Part I), 521-525 (Part II), 1397-1412 (Part III).

Sen, P.K. (1968). Estimates of Regression Coefficient Based on Kendall's tau. J. Am. Stat. Ass. 63, 324, 1379-1389.

Siegel, A.F. (1982). Robust Regression Using Repeated Medians. Biometrika, 69, 1, 242-244.
}
\author{
Lukasz Komsta, some fixes by Sven Garbade
}
\note{
This function should have compatibility with all 'lm' methods, but it is not guaranteed that they will work or have any cognitive value (this method is nonparametric). The compatibility was only introduced to use some basic methods from 'lm' without programming new functions.
}

\examples{
set.seed(1234)
x <- 1:100+rnorm(100)
y <- x+rnorm(100)
y[100] <- 200
fit <- mblm(y~x)
fit
summary(fit)
fit2 <- lm(y~x)
plot(x,y)
abline(fit)
abline(fit2,lty=2)
plot(fit)
residuals(fit)
fitted(fit)
plot(density(fit$slopes))
plot(density(fit$intercepts))
anova(fit)
anova(fit2)
anova(fit,fit2)
confint(fit)
AIC(fit,fit2)
}
\keyword{mblm}