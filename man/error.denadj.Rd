\name{error.denadj}
\alias{error.denadj}
\alias{error.cdfadj}
\title{
Compute the probability density function and cumulative probability density function of error, using localised bandwidths of residuals
}
\description{
With the estimated bandwidths of residuals, error density can be approximated by the kernel density estimator.
}
\usage{
error.denadj(ban, badj, eps, res.data)
}
\arguments{
  \item{ban}{Bandwidth of residuals}
  \item{badj}{Bandwidth adjustment factor}
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
 \code{\link{error.den}}
}
\keyword{method}
