\name{imse}
\alias{imse}
\title{IMSE Criterion}
\description{
Integrated Mean Squared Error:
\deqn{\mathrm{IMSE}(D;\theta)=\int_0^1(1-r(x)'R^{-1}r(x))\,dx}
}
\usage{
imse(D, theta, lambda = 1e-6)
}
\arguments{
  \item{D}{Design points (numeric vector).}
  \item{theta}{Correlation length parameter.}
  \item{lambda}{Regularization. Default \code{1e-6}.}
}
\value{Scalar IMSE value.}
