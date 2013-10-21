\name{specurves}
\alias{specurves}
\alias{fat}
\alias{protein}
\alias{moisture}
\docType{data}
\title{
Spectroscopy tecator data
}
\description{
This data set is a part of the original one which can be found at http://lib.stat/cmu.edu/datasets/tecator.
}
\usage{data(specurves)}
\format{
specurves: 215 by 100 matrix
fat: a vector of length 215
protein: a vector of length 215
moisture: a vector of length 215
}
\details{
For each unit, we observe one spectrometric curve which corresponds to the absorbance measured at 100 wavelengths (from 852 to 1050 in step of 2nm). For each measurement, we have at hand its fat content obtained by an analytic chemical processing
}
\source{
Nonparametric Functional Data Analysis website at http://www.lsp.ups-tlse.fr/staph/npfda/
}
\references{
C. Goutis (1998) "Second-derivative functional regression with applications to near infra-red spectroscopy",

F. Ferraty and P. Vieu (2002) "The functional nonparametric model and application to spectrometric

F. Ferraty and P. Vieu (2003) "Curve discrimination: A nonparametric functional approach", Computational

F. Ferraty and P. Vieu (2003) "Functional nonparametric statistics: A double infinite dimensional

F. Rossi and N. Delannay and B. Conan-Guez and M. Verleysen (2005) "Representation of functional

F. Ferraty and P. Vieu (2007) Nonparametric functional data analysis, New York: Springer.

H. Matsui and Y. Araki and S. Konishi (2008) "Multivariate regression modeling for functional
}
\examples{
data(specurves)
}
\keyword{datasets}