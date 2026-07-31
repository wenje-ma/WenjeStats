\name{sqrt_fit}
\alias{sqrt_fit}
\title{Fit Square-Root Posterior Approximation}
\description{
Fit the square-root basis expansion for an unnormalized posterior
density (Joseph, 2013). Given design points and evaluations, computes
the Gaussian correlation matrix, solves for basis coefficients, and
builds the pairwise matrices needed for marginal density computation.
}
\usage{
sqrt_fit(nu, hval, s = NULL, lambda = 1e-4, ...)
}
\arguments{
  \item{nu}{A numeric \eqn{n \times p} matrix of design points.}
  \item{hval}{A numeric vector of length \eqn{n} with the unnormalized
    posterior evaluated at each design point.}
  \item{s}{Optional length scales. If \code{NULL}, estimated automatically.}
  \item{lambda}{Regularization constant. Default \code{1e-4}.}
  \item{...}{Arguments passed to \code{estimate_scale} when \code{s = NULL}.}
}
\value{
An object of S3 class \code{"sqpost"} with components \code{nu},
\code{hval}, \code{s}, \code{coef}, \code{G}, \code{d}, and \code{M}.
}
\references{
Joseph, V. R. (2013). "Bayesian computation using design of experiments
and interpolation." \emph{Technometrics}, 55(4), 407-418.
}
\examples{
\dontrun{
fit <- sqrt_fit(nu, hev)
}
}
