\name{funopare.kernel}
\alias{funopare.kernel}
\title{
Functional Nadaraya-Watson estimator 
}
\description{
It implements the functional Nadaraya-Watson estimator to estimate the regression function. It depends on the type of semi-metric used as well as the optimal selection of bandwidth parameter
}
\usage{
funopare.kernel(Response, CURVES, PRED, bandwidth, ..., kind.of.kernel = "quadratic", 
	semimetric = "deriv")
}
\arguments{
  \item{Response}{A real-valued scalar response of length n}
  \item{CURVES}{An (n by p) matrix of discretised data of functional curves}
  \item{PRED}{An (n by k) matrix of discretised data of functional curves. \code{PRED} can be the same as the \code{CURVES} or the discretised data points of a new functional curve}
  \item{bandwidth}{A real-valued bandwidth parameter}
  \item{\dots}{Other arguments}
  \item{kind.of.kernel}{Type of kernel function. By default, it is the Epanechnikov kernel}
  \item{semimetric}{Type of semi-metric. By default, it is the semi-metric based on the qth order derivative, where q is an integer}
}
\details{
The functional NW estimator of the conditional mean can be expressed as a weighted average of response variable:
\eqn{\sum^n_{i=1}K_h(d(x_i,x))y_i/\sum^n_{i=1}K_h(d(x_i,x))},
where \eqn{K(\cdot)} is a kernel function which integrates to one, it has continuous derivative on the function support range. The semi-metric \eqn{d} is used to measure distances among curves. For a set of smooth curves, the semi-metric based on derivative should be considered. For a set of rough curves, the semi-metric based on functional principal components should be used. The bandwidth \eqn{h} controls the tradeoff between squared bias and variance in the mean squared error
}
\value{
\item{NWweit}{Estimated Nadaraya-Watson weights}
\item{Estimated.values}{Estimated values of the regression function}
\item{Predicted.values}{Predicted values of the regression function}
\item{band}{Bandwidth of the functional NW estimator}
\item{Mse}{In-sample mean squared error}
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, Computational Statistics, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, Computational Statistics and Data Analysis, 67, 185-198.

X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-
price density estimation, Journal of Econometrics, 153, 21-32.

F. Ferraty, I. Van Keilegom  and P. Vieu (2010) On the validity of the bootstrap in non-parametric functional regression, Scandinavian Journal of Statistics, 37, 286-306.

F. Ferraty and P. Vieu (2006) Nonparametric Functional Data Analysis: Theory and Practice, Springer, New York.

F. Ferraty and P. Vieu (2002) The functional nonparametric model and application to spectrometric data, Computational Statistics, 17, 545-564.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{bayMCMC_np_global}}, \code{\link{bayMCMC_np_local}}, \code{\link{bayMCMC_semi_global}}, \code{\link{bayMCMC_semi_local}}
}
\examples{
funopare.kernel(Response = simresp_np_normerr, CURVES = simcurve_smooth_normerr, 
	PRED = simcurve_smooth_normerr, bandwidth = 2.0, range.grid=c(0,pi), q=2, nknot=20)
}
\keyword{method}
