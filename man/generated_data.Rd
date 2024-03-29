% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{generated_data}
\alias{generated_data}
\title{Artificial genetic & phenotype data, for purposes of package tests}
\format{
A data frame with 10000 rows and 1 variables:
\describe{
\item{g}{gene count of 0, 1 or 2. Distibuted binomial with n=2 and p =0.3}
\item{u}{"unmeasured" confounder N(0,1)}
\item{errorX}{Error term for X,  ~exp(1)}
\item{errorY}{Error term for Y, ~N(0.1)}
\item{X}{Exposure. X= 2+ 0.25*g + u + errorX}
\item{linear.Y}{Linear outcome. Y =  X + 0.8u + errorY}
\item{quadratic.Y}{Quadratic outcome. Y= 2X^2 X + 0.8u + error Y}
\item{sqrt.Y}{Square root outcome. \eqn{Y = \sqrt{X} + 0.8*u+ errorY}}
\item{log.Y}{Log outcome. \eqn{Y = \log(X) + 0.8*U +errorY}}
\item{threshold.Y}{\eqn{X + 0.8* U +errorY} if \eqn{X>2}
and \eqn{0.8U + errorY} otherwise}
\item{fracpoly.Y}{Fractional polynomial \eqn{Y = X + 2/X + 0.8U + errorY}}
}
}
\usage{
generated_data
}
\description{
Artificial genetic & phenotype data, for purposes of package tests
}
\keyword{datasets}
