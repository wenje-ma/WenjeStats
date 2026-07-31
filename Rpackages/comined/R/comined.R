#' Logarithmic Distance (S-norm)
#'
#' Compute the logarithm of the S-norm of a vector, or the sum of logs
#' when \eqn{s = 0} (product distance).
#'
#' @param x A numeric vector.
#' @param s Distance exponent parameter.
#'   \itemize{
#'     \item If \eqn{s > 0}: returns \eqn{\frac{p}{s}\log\sum_{i=1}^p |x_i|^s},
#'           a scaled log of the S-norm.
#'     \item If \eqn{s = 0}: returns \eqn{\sum_{i=1}^p \log|x_i|},
#'           the log of the product of absolute values.
#'   }
#'
#' @return A numeric value: the log-distance.
#' @keywords internal
#' @examples
#' logdis(c(0.1, 0.2), s = 2)
#' logdis(c(0.1, 0.2), s = 0)
logdis <- function(x, s = 2) {
  if (s > 0) {
    length(x) / s * log(sum(abs(x)^s))
  } else {
    sum(log(abs(x)))
  }
}


#' Constrained Minimum Energy Design (COMiNED)
#'
#' Generate candidate points within a constrained region by combining
#' probit relaxation of inequality constraints with minimum energy
#' design. The algorithm uses a sequence of increasing penalty
#' parameters \eqn{\tau_0 < \tau_1 < \dots < \tau_T} to progressively
#' tighten the constraints.
#'
#' @section Mathematical Background:
#'
#' The feasible region is defined by \eqn{K} inequality constraints:
#' \deqn{\mathcal{X} = \{\boldsymbol{x} \in [0,1]^p : g_k(\boldsymbol{x}) \le 0,\; k = 1,\dots,K\}}
#'
#' Constraints are converted to a probability density via probit relaxation
#' (Golchi & Loeppky, 2015):
#' \deqn{f_\tau(\boldsymbol{x}) \propto \prod_{k=1}^K \Phi(-\tau g_k(\boldsymbol{x}))}
#'
#' where \eqn{\Phi} is the standard normal CDF. As \eqn{\tau \to \infty},
#' \eqn{f_\tau} approaches the indicator of the feasible region.
#'
#' At each stage \eqn{t}, the algorithm solves:
#' \deqn{\boldsymbol{x}_{m+1}^t = \arg\max_{\boldsymbol{x}} \min_{i=1:m}
#'   \left( \sum_{k=1}^K \log\Phi(-\tau_t g_k(\boldsymbol{x}))
#'         + \sum_{k=1}^K \log\Phi(-\tau_t g_k(\boldsymbol{x}_i))
#'         + 2p \log\|\boldsymbol{x}_i - \boldsymbol{x}\|_\alpha \right)}
#'
#' See Huang et al. (2021) for full details.
#'
#' @param n Number of design points to select.
#' @param p Number of input dimensions.
#' @param tau Numeric vector of strictly increasing penalty parameters.
#'   Typically starts at 0 and ends with a large value (e.g., \code{1e6}).
#' @param constraint A function that takes a numeric vector \code{x}
#'   (length \code{p}) and returns a numeric vector of constraint
#'   evaluations \eqn{g_1(\boldsymbol{x}), \dots, g_K(\boldsymbol{x})}.
#'   A point is feasible iff all elements are \eqn{\le 0}.
#' @param n.aug Integer. Augmentation multiplier; the initial lattice
#'   has \code{n * n.aug} points.
#' @param auto.scale Logical. If \code{TRUE}, automatically scale each
#'   constraint by its median absolute deviation (MAD) to balance
#'   constraints with different magnitudes.
#' @param s Distance exponent for the minimum energy design criterion.
#'   \itemize{
#'     \item \code{s = 2} (default): Euclidean distance.
#'     \item \code{s = 0}: Product (MaxPro-style) distance.
#'   }
#' @param gamma Numeric parameter passed to \code{\link[mined]{SelectMinED}}
#'   controlling the proportion of boundary points retained. Default is 0.5.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{med}}{Matrix (\code{n} x \code{p}) of final selected
#'       minimum energy design points.}
#'     \item{\code{med.lf}}{Numeric vector of log-density values at the
#'       selected points.}
#'     \item{\code{cand}}{Matrix of all candidate points generated across
#'       all stages.}
#'     \item{\code{cand.lf}}{Numeric vector of log-density values for all
#'       candidate points (at the final \eqn{\tau}).}
#'     \item{\code{cand.min.dist}}{Minimum Euclidean distance among
#'       candidate points.}
#'     \item{\code{cand.min.ddist}}{Minimum one-dimensional projection
#'       distance among candidate points.}
#'     \item{\code{feasible.idx}}{Logical vector indicating which rows of
#'       \code{cand} satisfy all constraints (\eqn{g_k \le 0}).}
#'   }
#'
#' @references
#' \itemize{
#'   \item Huang, C., Joseph, V. R., & Mak, S. (2021). "Constrained
#'         minimum energy designs." \emph{Technometrics}, 63(4), 445-458.
#'   \item Golchi, S. & Loeppky, J. L. (2015). "A probabilistic
#'         relaxation approach for constrained space-filling designs."
#'         \emph{arXiv preprint} arXiv:1508.06861.
#'   \item Joseph, V. R., Dasgupta, T., Tuo, R., & Wu, C. F. J. (2015).
#'         "Sequential minimum energy design." \emph{Journal of the
#'         American Statistical Association}, 110(511), 1116-1130.
#' }
#'
#' @importFrom mined Lattice
#' @importFrom mined SelectMinED
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define constraint: g1, g2, g3 from Equation (4.30)
#' constraint <- function(x) {
#'   c1 <- x[1] - sqrt(50 * (x[2] - 0.52)^2 + 2) + 1
#'   c2 <- sqrt(120 * ((x[2] - 0.48)^2 + 1)) - 0.75 - x[1]
#'   c3 <- 0.65^2 - x[1]^2 - x[2]^2
#'   c(c1, c2, c3)
#' }
#'
#' tau <- c(0, exp(1:7), 1e6)
#' result <- comined(n = 53, p = 2, tau = tau,
#'                   constraint = constraint, n.aug = 5)
#'
#' # Extract feasible points
#' feasible <- result$cand[result$feasible.idx, ]
#' cat("Feasible points:", sum(result$feasible.idx), "/", nrow(result$cand), "\n")
#' }
comined <- function(n, p, tau, constraint, n.aug,
                    auto.scale = FALSE, s = 2, gamma = 0.5) {

  # ---- Stage 1: Initial lattice design ----
  samp <- Lattice(n * n.aug, p)
  min.dist <- min(dist(samp))
  min.ddist <- Inf
  for (i in 1:p) {
    ddist <- min(dist(samp[, i]))
    if (ddist < min.ddist) min.ddist <- ddist
  }

  # ---- Evaluate constraints ----
  samp.gval <- matrix(t(apply(samp, 1, constraint)), nrow = nrow(samp))

  # ---- Stage 2: First MED selection (tau[2]) ----
  scale <- rep(1, ncol(samp.gval))
  if (auto.scale) scale <- apply(samp.gval, 2, mad, center = 0)

  samp.lf <- apply(samp.gval, 1, function(x) {
    sum(pnorm(-tau[2] * x / scale, log.p = TRUE))
  })

  med.op <- SelectMinED(samp, samp.lf, n, gamma = gamma, s = s)
  samp.med <- med.op$points

  # ---- Stages 3..T: progressively tighten constraints ----
  for (k in 3:length(tau)) {
    min.dist <- min.dist / 2
    min.ddist <- min.ddist / 2

    # Determine rounding precision from min.ddist
    no.decimal <- attr(
      regexpr("(?<=\\.)0+", format(min.ddist, scientific = FALSE), perl = TRUE),
      "match.length"
    ) + 1

    # Adaptive augmentation around current MED points
    samp.med.dist <- as.matrix(dist(samp.med))
    samp.aug <- NULL
    for (i in 1:n) {
      nn.idx <- order(samp.med.dist[i, ])[2:(n.aug + 1)]
      # Midpoints
      samp.aug <- rbind(
        samp.aug,
        0.5 * (samp.med[nn.idx, ] + rep(1, n.aug) %*% t(samp.med[i, ]))
      )
      # Extended points
      samp.aug <- rbind(
        samp.aug,
        0.5 * (3 * samp.med[nn.idx, ] - rep(1, n.aug) %*% t(samp.med[i, ]))
      )
    }

    # Remove duplicates
    samp.aug.rep <- round(samp.aug, digits = no.decimal)
    samp.aug <- samp.aug[!duplicated(samp.aug.rep), ]

    # Remove points outside [0,1]^p
    samp.aug.out <- apply(samp.aug, 1, function(x) any(x < 0 | x > 1))
    samp.aug <- samp.aug[!samp.aug.out, ]

    # Remove points already in samp
    no.aug <- nrow(samp.aug)
    if (no.aug > 0) {
      samp.rep <- rbind(samp.aug, samp)
      samp.rep <- round(samp.rep, digits = no.decimal)
      samp.aug <- samp.aug[!duplicated(samp.rep, fromLast = TRUE)[1:no.aug], ]
    }

    # Evaluate constraints on augmented points and merge
    if (nrow(samp.aug) > 0) {
      samp.aug.gval <- matrix(
        t(apply(samp.aug, 1, constraint)), nrow = nrow(samp.aug)
      )
      samp <- rbind(samp, samp.aug)
      samp.gval <- rbind(samp.gval, samp.aug.gval)
    }

    # Update scale if auto-scaling
    if (auto.scale) scale <- apply(samp.gval, 2, mad, center = 0)

    # Compute log-density with current tau[k]
    samp.lf <- apply(samp.gval, 1, function(x) {
      sum(pnorm(-tau[k] * x / scale, log.p = TRUE))
    })

    # Estimate minimum log-distance for threshold
    if (s == 0) {
      hp <- floor(p / 2)
      nl <- sqrt((min.dist^2 - min.ddist^2 * hp) / (p - hp))
      min.logdis <- logdis(c(rep(min.ddist, hp), rep(nl, p - hp)), s = s)
    } else {
      min.logdis <- logdis(rep(min.dist / sqrt(p), p), s = s)
    }

    # Threshold: keep points with high log-density
    samp.lf.cv <- sort(samp.lf, decreasing = TRUE)[n] + 2.5 * min.logdis
    samp.cand <- samp[samp.lf > samp.lf.cv, , drop = FALSE]
    samp.cand.lf <- samp.lf[samp.lf > samp.lf.cv]

    # Select n MED points from the filtered candidates
    med.op <- SelectMinED(samp.cand, samp.cand.lf, n,
                          gamma = gamma, s = s)
    samp.med <- med.op$points
  }

  # ---- Final results ----
  samp.med.lf <- med.op$logf
  feasible.idx <- !apply(samp.gval, 1, function(x) any(x > 0))

  list(
    med           = samp.med,
    med.lf        = samp.med.lf,
    cand          = samp,
    cand.lf       = samp.lf,
    cand.min.dist = min.dist,
    cand.min.ddist = min.ddist,
    feasible.idx  = feasible.idx
  )
}
