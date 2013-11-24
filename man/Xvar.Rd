\name{Xvar}
\alias{Xvar}
\alias{tau_semierr}
\alias{simresp_semi_normerr}
\docType{data}
\title{
Simulated real-valued predictors in the semi-functional partial linear model
}
\description{
Simulated real-valued predictors in the semi-functional partial linear model. It is a 50 by 2 matrix, where the column variables are generated from \eqn{U[0,1]} and the true regression coefficients are (-1,2). Note that the estimation of the regression coefficient for these predictors depends crucially on the bandwidth parameter estimated in the functional Nadaraya-Watson estimator of the regression function
}
\usage{data(Xvar)}
\source{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.
}
\examples{
data(Xvar)
data(tau_semierr)
data(simresp_semi_normerr)
}
\keyword{datasets}
