\name{logpriors_admkr}
\alias{logpriors_admkr}
\title{
Compute the marginal likelihood using Chib's (1995) method
}
\description{
The log marginal likelihood can be computed as the log likelihood + log prior - log posterior. The latter one can be estimated from the MCMC.
}
\usage{
logpriors_admkr(h2, prior_alpha, prior_beta)
}
\arguments{
  \item{h2}{Squared bandwidths}
  \item{prior_alpha}{Hyper-parameter of the inverse gamma distribution}
  \item{prior_beta}{Hyper-parameter of the inverse gamma distribution}
}
\value{
The value of log prior density in the marginal likelihood
}
\references{
S. Chib (1995) `Marginal likelihood from the Gibbs output', \emph{Journal of the American Statistical Association}, \bold{90}(432), 1313-1321.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{logdensity_admkr}}, \code{\link{loglikelihood_global_admkr}}
}
\keyword{method}
