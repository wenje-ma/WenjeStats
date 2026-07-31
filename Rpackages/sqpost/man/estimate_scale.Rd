\name{estimate_scale}
\alias{estimate_scale}
\title{Estimate Length Scales via Cross-Validation}
\description{
Estimate the length scale parameters \eqn{s_1,\dots,s_p} for the
square-root posterior approximation by minimizing the log of the
weighted mean squared cross-validation error (wMSCV).
}
\usage{
estimate_scale(nu, hval, lambda = 1e-4,
               lower = NULL, upper = NULL, ...)
}
\arguments{
  \item{nu}{A numeric \eqn{n \times p} matrix of design points.}
  \item{hval}{A numeric vector of length \eqn{n} giving the unnormalized
    posterior density evaluated at each design point.}
  \item{lambda}{Regularization constant. Default \code{1e-4}.}
  \item{lower, upper}{Bounds for the optimization. Default to
    \code{ini.s/100} and \code{100*ini.s}.}
  \item{...}{Additional arguments passed to \code{optim}.}
}
\value{
A numeric vector of length \eqn{p} of estimated length scales.
}
\references{
Joseph, V. R. (2013). "Bayesian computation using design of experiments
and interpolation." \emph{Technometrics}, 55(4), 407-418.
}
\examples{
\dontrun{
s <- estimate_scale(nu, hev)
}
}
