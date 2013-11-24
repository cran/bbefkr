\name{error.den}
\alias{error.den}
\alias{admkr.cdf}
\title{
Compute the probability density function and cumulative probability density function of the error, using a global bandwidth of residuals 		
}
\description{
With the estimated bandwidth of residuals, error density can be approximated by the kernel density estimator.
}
\usage{
error.den(band, eps, res.data)
}
\arguments{
  \item{band}{Bandwidth of residuals}
  \item{eps}{Grid point}
  \item{res.data}{Residuals obtained from the estimated conditional mean}
}
\value{
Numerical values
}
\references{
H. L. Shang (2013) `Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density', \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

H. L. Shang (2013) `Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density', \emph{Computational Statistics}, in press.
}
\author{
Han Lin Shang
}
\seealso{
 \code{\link{error.denadj}}
}
\keyword{method}
