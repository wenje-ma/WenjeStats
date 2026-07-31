\name{gauss_corr}
\alias{gauss_corr}
\title{Gaussian Correlation Function (1D)}
\description{
Compute the Gaussian correlation between two sets of 1D points:
\deqn{R(x_i, y_j) = \exp\left\{-\left(\frac{x_i - y_j}{\theta}\right)^2\right\}}
}
\usage{
gauss_corr(x, y = NULL, theta)
}
\arguments{
  \item{x}{Numeric vector of length \eqn{n}.}
  \item{y}{Numeric vector of length \eqn{m}. If \code{NULL}, uses \code{x}.}
  \item{theta}{Positive numeric correlation length parameter.}
}
\value{An \eqn{n \times m} matrix.}
\examples{
gauss_corr(0:5/5, theta = 0.2)
}
