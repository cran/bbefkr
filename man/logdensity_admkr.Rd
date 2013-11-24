\name{logdensity_admkr}
\alias{logdensity_admkr}
\title{
Compute the marginal likelihood using Chib's (1995) method
}
\description{
The log marginal likelihood can be computed as the log likelihood + log prior - log posterior. The latter one can be estimated from the MCMC.
}
\usage{
logdensity_admkr(tau2, cpost)
}
\arguments{
  \item{tau2}{Squared bandwidths recorded at each iteration}
  \item{cpost}{Squared bandwidths recorded at only the iteration defined by the variable \code{step}}
}
\value{
The value of the estimated posterior density
}
\references{
S. Chib (1995) `Marginal likelihood from the Gibbs output', \emph{Journal of the American Statistical Association}, \bold{90}(432), 1313-1321.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{logpriors_admkr}}, \code{\link{loglikelihood_admkr}}
}
\keyword{method}
