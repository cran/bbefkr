\name{simulate_error}
\alias{simulate_error}
\title{
Simulate errors
}
\description{
Simulate errors which has the specified error density
}
\usage{
simulate_error(samplesize, errordensity)
}
\arguments{
  \item{samplesize}{Sample size}
  \item{errordensity}{Type of error density}
}
\value{
A vector of simulated error, which follows the specific error density
}
\references{
J. S. Marron and M. P. Wand (1992) `Exact Mean Integrated Squared Error', Annals of Statistics, 20(2), 712-736.
}
\author{
Han Lin Shang
}
\examples{
simulate_error(samplesize = 100, errordensity = "normal")
}
\keyword{method}
