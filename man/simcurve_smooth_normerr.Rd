\name{simcurve_smooth_normerr}
\alias{simcurve_smooth_normerr}
\alias{simcurve_rough_normerr}
\alias{simresp_np_normerr}
\alias{tau_normerr}
\docType{data}
\title{Simulated data set}
\description{
Simulated data
}
\usage{data(simcurve_smooth_normerr)}
\format{
\code{simcurve_smooth_normerr}: 50 by 100

\code{simcurve_rough_normerr}: 50 by 100

\code{simresp_normerr}: 50 by 1

\code{tau_normerr}: 50 by 1 	
}
\details{
	The simulated discretised curves are defined as
	\eqn{x_i(t_j) = a_i cos(2t_j)+b_isin(4t_j)+c_i(t_j^2-\pi \times t_j+2/9\pi^2)},
	where t represents the function support range and \eqn{0\leq t_1\leq t_2\dots\leq \pi} are equispaced points within the function support range, \eqn{a_i}, \eqn{b_i} and \eqn{c_i} are independently drawn from a uniform distribution on [0,1], and \eqn{n} represents the sample size. For simulating a set of rough curves, we add one extra term \eqn{d_j} generated from \eqn{U(-0.1, 0.1)}.
	Having defined functional curves, we then construct the regression mean function \eqn{\tau=10\times (a_i^2-b_i^2)}. Then the response variable is obtained by adding the regression mean function with a set of errors generated from a standard normal distribution
}
\source{
H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

F. Ferraty, I. Van Keilegom, P. Vieu (2010) On the validity of the bootstrap in non-parametric functional regression, \emph{Scandinavian Journal of Statistics}, \bold{37}(2), 286-306.
}
\examples{
data(simcurve_normerr)
}
\keyword{datasets}
