\name{SIF}
\alias{SIF}
\title{
Simulation inefficiency factor
}
\description{
Calculate the simulation inefficiency factor, as a means of checking the convergence of MCMC
}
\usage{
SIF(BAND_MATRIX, NUM_ITERATIONS, NUM_BATCH)
}
\arguments{
  \item{BAND_MATRIX}{Bandwidth parameters drawn from the MCMC}
  \item{NUM_ITERATIONS}{Number of total iterations}
  \item{NUM_BATCH}{Number of batch samples}
}
\value{
Numerical values
}
\references{
S. Kim, N. Shephard and S. Chib (1998) `Stochastic volatility: likelihood inference and comparison with ARCH models', \emph{The Review of Economic Studies}, \bold{65}(3), 361-393.

R. Meyer and J. Yu (2000) `BUGS for a Bayesian analysis of stochastic volatility models', \emph{Econometrics Journal}, \bold{3}(2), 198-215.

Y. K. Tse, X. Zhang and J. Yu (2004) `Estimation of hyperbolic diffusion using the Markov chain Monte Carlo method', \emph{Quantitative finance}, \bold{4}(2), 158-169.

X. Zhang, R. D. Brooks and M. L. King (2009) `A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation', \emph{Journal of Econometrics}, \bold{153}(1), 21-32.
}
\author{
Han Lin Shang
}
\keyword{method}
