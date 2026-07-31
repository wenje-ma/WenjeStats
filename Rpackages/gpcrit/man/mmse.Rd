\name{mmse}
\alias{mmse}
\title{MMSE Criterion}
\description{
Maximum Mean Squared Error:
\deqn{\mathrm{MMSE}(D;\theta)=\max_{x\in[0,1]}(1-r(x)'R^{-1}r(x))}
}
\usage{
mmse(D, theta, n_grid = 301, lambda = 1e-6)
}
\arguments{
  \item{D}{Design points (numeric vector).}
  \item{theta}{Correlation length parameter.}
  \item{n_grid}{Grid size for max approximation. Default 301.}
  \item{lambda}{Regularization. Default \code{1e-6}.}
}
\value{Scalar MMSE value.}
