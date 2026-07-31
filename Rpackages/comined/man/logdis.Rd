\name{logdis}
\alias{logdis}
\title{Logarithmic Distance (S-norm)}
\description{
Compute the logarithm of the S-norm of a vector, or the sum of logs
when \eqn{s = 0} (product distance).
}
\usage{
logdis(x, s = 2)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{s}{Distance exponent parameter.
    \itemize{
      \item If \eqn{s > 0}: returns \eqn{\frac{p}{s}\log\sum_{i=1}^p |x_i|^s}.
      \item If \eqn{s = 0}: returns \eqn{\sum_{i=1}^p \log|x_i|}.
    }
  }
}
\value{
A numeric value: the log-distance.
}
\keyword{internal}
