\name{marginal}
\alias{marginal}
\title{Marginal Posterior Density}
\description{
Compute the marginal posterior density for one dimension from a
fitted \code{sqpost} object using the closed-form expression
(Joseph, 2013, Eq. 4.53).
}
\usage{
marginal(object, grid, dim = 1)
}
\arguments{
  \item{object}{An object of class \code{"sqpost"} from \code{sqrt_fit}.}
  \item{grid}{A numeric vector of grid points for evaluation.}
  \item{dim}{Integer. Which dimension (column index) to marginalize.
    Default is 1.}
}
\value{
A numeric vector of marginal density values at the grid points.
}
\references{
Joseph, V. R. (2013). "Bayesian computation using design of experiments
and interpolation." \emph{Technometrics}, 55(4), 407-418.
}
\examples{
\dontrun{
fit <- sqrt_fit(nu, hev)
den1 <- marginal(fit, seq(0, 1, length.out = 100), dim = 1)
}
}
