\name{gauss_corr}
\alias{gauss_corr}
\title{Gaussian Correlation Matrix}
\description{
Build the Gaussian correlation matrix for a set of points with
dimension-specific length scales.
}
\usage{
gauss_corr(x, s)
}
\arguments{
  \item{x}{A numeric \eqn{n \times p} matrix.}
  \item{s}{A numeric vector of length \eqn{p} of positive length scales.}
}
\value{
An \eqn{n \times n} symmetric correlation matrix.
}
\keyword{internal}
