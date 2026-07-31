\name{gp_rmse}
\alias{gp_rmse}
\title{GP RMSE}
\description{
Compute standardized RMSE at test points given a design and correlation
parameter. Formula: \eqn{s(x)=\sqrt{1-r(x)'R^{-1}r(x)}}.
}
\usage{
gp_rmse(D, theta, test, lambda = 1e-6)
}
\arguments{
  \item{D}{Design points (numeric vector).}
  \item{theta}{Correlation length parameter.}
  \item{test}{Test points (numeric vector).}
  \item{lambda}{Regularization for R inverse. Default \code{1e-6}.}
}
\value{Numeric vector of RMSE values.}
