\name{comined}
\alias{comined}
\title{Constrained Minimum Energy Design (COMiNED)}
\description{
Generate candidate points within a constrained region by combining
probit relaxation of inequality constraints with minimum energy
design. The algorithm uses a sequence of increasing penalty
parameters \eqn{\tau_0 < \tau_1 < \dots < \tau_T} to progressively
tighten the constraints.
}
\usage{
comined(n, p, tau, constraint, n.aug,
        auto.scale = FALSE, s = 2, gamma = 0.5)
}
\arguments{
  \item{n}{Number of design points to select.}
  \item{p}{Number of input dimensions.}
  \item{tau}{Numeric vector of strictly increasing penalty parameters.
    Typically starts at 0 and ends with a large value (e.g., \code{1e6}).}
  \item{constraint}{A function that takes a numeric vector \code{x}
    (length \code{p}) and returns a numeric vector of constraint
    evaluations. A point is feasible iff all elements are \eqn{\le 0}.}
  \item{n.aug}{Integer. Augmentation multiplier; the initial lattice
    has \code{n * n.aug} points.}
  \item{auto.scale}{Logical. If \code{TRUE}, automatically scale each
    constraint by its median absolute deviation (MAD).}
  \item{s}{Distance exponent for the minimum energy design criterion.
    \code{s = 2} (default): Euclidean distance; \code{s = 0}: Product
    (MaxPro-style) distance.}
  \item{gamma}{Numeric parameter controlling the proportion of boundary
    points retained in \code{\link[mined]{SelectMinED}}. Default 0.5.}
}
\value{
A list with components:
  \describe{
    \item{\code{med}}{Matrix (\code{n} x \code{p}) of final selected MED points.}
    \item{\code{med.lf}}{Numeric vector of log-density values at the selected points.}
    \item{\code{cand}}{Matrix of all candidate points generated across all stages.}
    \item{\code{cand.lf}}{Numeric vector of log-density values for all candidates.}
    \item{\code{cand.min.dist}}{Minimum Euclidean distance among candidates.}
    \item{\code{cand.min.ddist}}{Minimum one-dimensional projection distance.}
    \item{\code{feasible.idx}}{Logical vector: which rows of \code{cand} satisfy all constraints.}
  }
}
\references{
Huang, C., Joseph, V. R., & Mak, S. (2021). "Constrained minimum energy designs."
\emph{Technometrics}, 63(4), 445-458.

Golchi, S. & Loeppky, J. L. (2015). "A probabilistic relaxation approach for
constrained space-filling designs." \emph{arXiv preprint} arXiv:1508.06861.

Joseph, V. R., Dasgupta, T., Tuo, R., & Wu, C. F. J. (2015). "Sequential minimum
energy design." \emph{Journal of the American Statistical Association}, 110(511), 1116-1130.
}
\examples{
\dontrun{
constraint <- function(x) {
  c1 <- x[1] - sqrt(50 * (x[2] - 0.52)^2 + 2) + 1
  c2 <- sqrt(120 * ((x[2] - 0.48)^2 + 1)) - 0.75 - x[1]
  c3 <- 0.65^2 - x[1]^2 - x[2]^2
  c(c1, c2, c3)
}
tau <- c(0, exp(1:7), 1e6)
result <- comined(53, 2, tau, constraint, n.aug = 5)
feasible <- result$cand[result$feasible.idx, ]
}
}
