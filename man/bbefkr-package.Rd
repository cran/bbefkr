\name{bbefkr-package}
\alias{bbefkr-package}
\alias{bbefkr}
\docType{package}
\title{
Bayesian bandwidth estimation for the functional kernel regression with unknown error density
}
\description{
This package aims to estimate bandwidths in the regression function and kernel-form error density simultaneously, using a Bayesian sampling algorithm. We demonstrate this Bayesian sampling algorithm using a functional nonparametric regression and a semi-functional partial linear model
}
\details{
The regression function is approximated by the functional Nadaraya-Watson estimator, while the unknown error density is approximated by a kernel density of residuals. In both regression function and error density, they depend crucially on the selection of the optimal bandwidths.
}
\author{
Han Lin Shang

Maintainer: Han Lin Shang <H.Shang@soton.ac.uk>
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, \emph{Journal of Econometrics}, \bold{153}, 21-32.
}
\keyword{ package }
\seealso{
\code{\link[bbefkr]{bayMCMC_np_global}}, \code{\link[bbefkr]{bayMCMC_np_local}}, \code{\link[bbefkr]{bayMCMC_semi_global}}, \code{\link[bbefkr]{bayMCMC_semi_local}}
}
