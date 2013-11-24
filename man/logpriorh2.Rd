\name{logpriorh2}
\alias{logpriorh2}
\title{
Prior density of the squared bandwidth parameters
}
\description{
Inverse gamma distribution is used as the prior density of bandwidth parameters
}
\usage{
logpriorh2(h2, prior_alpha, prior_beta)
}
\arguments{
  \item{h2}{Squared bandwidths}
  \item{prior_alpha}{Hyper-parameter of the inverse gamma distribution}
  \item{prior_beta}{Hyper-parameter of the inverse gamma distribution}
}
\value{
Prior density
}
\references{
H. L. Shang (2013) `Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density', \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

H. L. Shang (2013) `Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density', \emph{Computational Statistics}, in press.
}
\author{
Han Lin Shang
}
\keyword{method}
