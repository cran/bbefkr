\name{bayMCMC_np_local}
\alias{bayMCMC_np_local}
\title{
Bayesian bandwidth estimation for a functional nonparametric regression with homoscedastic errors
}
\description{
Estimate the bandwidths in the regression function approximated by the functional Nadaraya-Watson estimator and kernel-form error density with localised bandwidths, in a functional nonparametric regression
}
\usage{
bayMCMC_np_local(data_x, data_y, data_xnew, warm = 1000, M = 1000, 
	mutprob = 0.44, errorprob = 0.44, epsilonprob = 0.44, mutsizp = 1, 
		errorsizp = 1, epsilonsizp = 1, prior_alpha = 1, prior_beta = 0.05, 
			err_int = c(-10, 10), err_ngrid = 10001, num_batch = 20, 
				step = 10, alpha = 0.95, ...)
}
\arguments{
  \item{data_x}{An (n by p) matrix of discretised data points of functional curves}
  \item{data_y}{A scalar-valued response of length n}
  \item{data_xnew}{A matrix of discretised data points of new functional curve(s)}
  \item{warm}{Number of iterations for the burn-in period}
  \item{M}{Number of iterations for the Markov chain Monte Carlo (MCMC)}
  \item{mutprob}{Optimal acceptance rate of the random-walk Metropolis algorithm for sampling the bandwidth in the regression function}
  \item{errorprob}{Optimal acceptance rate of the random-walk Metropolis algorithm for sampling the bandwidth in the kernel-form error density}
  \item{epsilonprob}{Optimal acceptance rate of the random-walk Metropolis algorithm for sampling the bandwidth adjust factor in the kernel-form error density}
  \item{mutsizp}{Initial step length of the random-walk Metropolis algorithm for sampling the bandwidth in the regression function. Its value will be updated at each iteration to achieve the optimal acceptance rate, given the MCMC converges to its target distribution}
  \item{errorsizp}{Initial step length of the random-walk Metropolis algorithm for sampling the bandwidth in the kernel-form error density. Its value will be updated at each iteration to achieve the optimal acceptance rate, given the MCMC converges to its target distribution}
  \item{epsilonsizp}{Initial step length of the random-walk Metropolis algorithm for sampling the bandwidth adjustment factor in the kernel-form error density. Its value will be updated at each iteration to achieve the optimal acceptance rate}
  \item{prior_alpha}{Hyperparameter of the inverse gamma prior distribution for the squared bandwidths}
  \item{prior_beta}{Hyperparameter of the inverse gamma prior distribution for the squared bandwidths}
  \item{err_int}{Range of the error-density grid for computing the probability density function and cumulative probability density function}
  \item{err_ngrid}{Number of the error-density grid points}
  \item{num_batch}{Number of batches to assess the convergence of the MCMC}
  \item{step}{Thinning parameter. For example, when \code{step=10}, it keeps every 10th iteration of the MCMC output}
  \item{alpha}{The nominal coverage probability of the prediction interval, customarily 95 percent}  
  \item{...}{Other arguments used to define semi-metric. For a set of smoothed functional data, the semi-metric based on derivative is suggested. For a set of rough functional data, the semi-metric based on the functional principal component analysis is suggested}
}
\details{
The Bayesian method estimates the bandwidths in the regression function and kernel-form error density. It performs better than the functional cross validation in terms of estimation accuracy, since the latter one does not utilise the information about the unknown error density. Furthermore, it can estimate error density more accurate than the Bayesian method with a global bandwidth
}
\value{
\item{xpfinalres}{Estimated bandwidths}
\item{mhat}{Estimated regression function}
\item{sif_value}{Simulation inefficiency factor}
\item{mlikeres}{Marginal likelihood calculated using the Chib's (1995) method}
\item{acceptnwMCMC}{Acceptance rate for sampling bandwidth in the regression function}
\item{accepterroMCMC}{Acceptance rate for sampling bandwidth in the kernel-form error density}
\item{acceptepsilonMCMC}{Acceptance rate for sampling bandwidth adjustment factor in the kernel-form error density}
\item{fore.den.mkr}{Estimated probability density function of the error}
\item{fore.cdf.mkr}{Estimated cumulative density function of the error}
\item{point forecast}{Predicted response}
\item{PI}{Prediction interval of response}
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a functional nonparametric regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, 67, 185-198.

X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, \emph{Journal of Econometrics}, \bold{153}, 21-32.

F. Ferraty, I. Van Keilegom  and P. Vieu (2010) On the validity of the bootstrap in non-parametric functional regression, \emph{Scandinavian Journal of Statistics}, \bold{37}, 286-306.

R. Meyer and J. Yu (2000) BUGS for a Bayesian analysis of stochastic volatility models, \emph{Econometircs Journal}, \bold{3}(2), 198-215.

S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}(432), 1313-1321.
}
\author{
Han Lin Shang
}
\note{
It can be time-consuming when the sample size is large, say above 250
}

\seealso{
\code{\link{bayMCMC_np_global}}, \code{\link{bayMCMC_semi_local}}
}
\keyword{method}
