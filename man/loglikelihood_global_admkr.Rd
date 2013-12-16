\name{loglikelihood_global_admkr}
\alias{loglikelihood_global_admkr}
\alias{loglikelihood_local_admkr}
\title{
Compute the marginal likelihood using Chib's (1995) method
}
\description{
The log marginal likelihood can be computed as the log likelihood + log prior - log posterior. The latter one can be estimated from the MCMC.
}
\usage{
loglikelihood_global_admkr(h, resid)
}
\arguments{
  \item{h}{Estimated bandwidth}
  \item{resid}{Estimated residuals}
}
\value{
The value of log likelihood
}
\references{
S. Chib (1995) `Marginal likelihood from the Gibbs output', \emph{Journal of the American Statistical Association}, \bold{90}(432), 1313-1321.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{logpriors_admkr}}, \code{\link{logdensity_admkr}}
}
\keyword{method}
