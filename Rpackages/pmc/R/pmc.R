#' Population (Quasi-)Monte Carlo
#'
#' Runs the Population Monte Carlo (PMC) algorithm (Cappé et al., 2004) or
#' the Population Quasi-Monte Carlo (PQMC) algorithm (Huang et al., 2022)
#' for sampling from a target distribution known up to a normalising constant.
#'
#' The algorithm iterates \code{steps} rounds. In each round:
#' \enumerate{
#'   \item \strong{Sampling}: draw \code{J} samples from each of the \code{K}
#'         proposal distributions (Gaussian centred at current centres).
#'   \item \strong{Weighting}: compute importance weights using a deterministic
#'         mixture of the proposals.
#'   \item \strong{Adaptation}: resample \code{K} new centres from the weighted
#'         pool, and optionally adapt the scale parameter \code{sigma}.
#' }
#'
#' @param logf function computing the log target density (up to a constant).
#' @param K integer; number of proposal centres.
#' @param J integer; number of samples per proposal per iteration.
#' @param steps integer; number of PMC iterations.
#' @param init numeric matrix (\code{K x p}) of initial centres.
#' @param sampling character; \code{"random"} for standard PMC,
#'        \code{"qmc"} for quasi-Monte Carlo sampling (PQMC).
#' @param resampling character; resampling method passed to \code{\link{resample}}.
#' @param output character; \code{"weighted samples"} returns all samples and
#'        weights, \code{"estimator"} returns the integrated estimator.
#' @param sigma initial/ fixed scale parameter. If \code{NULL} and
#'        \code{sigma.adapt = TRUE}, initialised as \code{min(dist(init))}.
#' @param sigma.adapt logical; whether to adapt \code{sigma} each iteration.
#' @param visualization logical; if \code{TRUE} and \code{p = 2}, draw
#'        contour plots of the target at each iteration.
#'
#' @return If \code{output = "weighted samples"}, a list with components:
#'   \item{samp.all}{all samples (matrix).}
#'   \item{samp.all.logwts}{log importance weights of all samples.}
#'   \item{center.all}{all proposal centres (initial + adapted).}
#'   \item{ess}{effective sample size at each iteration.}
#'   If \code{output = "estimator"}, a list with components:
#'   \item{m.std}{standard PMC estimate of \eqn{E_f[h]}.}
#'   \item{m.wts}{weighted PMC estimate.}
#'   \item{m.las}{last-adapted estimate.}
#'   \item{z.std, z.wts}{normalising constant estimates.}
#'
#' @importFrom stats rnorm dnorm qnorm
#' @importFrom randtoolbox sobol
#' @export
#'
#' @examples
#' \dontrun{
#' logmix <- function(x) {
#'   log(1/5 * sum(dmvnorm(x, mean = c(0.5, 0.5), sigma = diag(0.1, 2))))
#' }
#' init <- matrix(runif(20), ncol = 2)
#' result <- pmc(logmix, K = 10, J = 5, steps = 3, init = init)
#' }
pmc <- function(logf, K, J, steps, init,
                sampling = c("random", "qmc"),
                resampling = c("multinomial", "residual", "systematic",
                               "stratified", "sp"),
                output = c("weighted samples", "estimator"),
                sigma = NULL, sigma.adapt = TRUE, visualization = FALSE) {

  sampling <- match.arg(sampling)
  resampling <- match.arg(resampling)
  output <- match.arg(output)

  if (nrow(init) != K) stop("nrow(init) must equal K")
  center <- init

  # initialise sigma
  if (is.null(sigma)) {
    if (sigma.adapt) {
      sigma <- min(dist(center))
    } else {
      stop("sigma must be provided if sigma.adapt = FALSE")
    }
  }

  p <- ncol(init)
  if (visualization && p != 2) visualization <- FALSE
  if (visualization) {
    x1 <- x2 <- seq(0, 1, length.out = 101)
    x.grid <- expand.grid(x1, x2)
    z <- matrix(exp(apply(x.grid, 1, logf)), 101, 101)
    graphics::contour(x1, x2, z, drawlabels = FALSE, nlevels = 15,
                      main = sprintf("t = %d", 0))
    graphics::points(center, pch = 16, cex = 1, col = "red")
  }

  samp.all <- NULL
  samp.all.logwts <- NULL
  center.all <- center
  ess <- rep(NA, steps)

  for (t in 1:steps) {
    # --- Sampling ---
    samp <- NULL
    if (sampling == "random") {
      for (i in 1:K) {
        noise <- matrix(stats::rnorm(J * p, mean = 0, sd = sigma), ncol = p)
        samp <- rbind(samp, rep(1, J) %*% t(center[i, ]) + noise)
      }
    } else {  # qmc
      for (i in 1:K) {
        noise <- stats::qnorm(
          suppressWarnings(randtoolbox::sobol(J, p, scrambling = 1,
                                              seed = sample(1e6, 1))),
          mean = 0, sd = sigma)
        samp <- rbind(samp, rep(1, J) %*% t(center[i, ]) + noise)
      }
    }

    # --- Weighting ---
    samp.dm.logwts <- matrix(0, nrow = nrow(samp), ncol = K)
    for (i in 1:K) {
      for (j in 1:nrow(samp)) {
        dist <- samp[j, ] - center[i, ]
        samp.dm.logwts[j, i] <- sum(stats::dnorm(dist, mean = 0, sd = sigma, log = TRUE))
      }
    }
    samp.logq <- apply(samp.dm.logwts, 1, logaddexp) - log(K)
    samp.logf <- apply(samp, 1, logf)
    samp.logf[is.na(samp.logf)] <- -Inf
    samp.logwts <- samp.logf - samp.logq
    samp.wts <- exp(samp.logwts - logaddexp(samp.logwts))

    # --- Store ---
    samp.all <- rbind(samp.all, samp)
    samp.all.logwts <- c(samp.all.logwts, samp.logwts)
    ess[t] <- 1 / sum(samp.wts^2)

    # --- Covariance adaptation ---
    if (sigma.adapt && t < steps) {
      dm.wts <- exp(samp.dm.logwts - (samp.logq + log(K)) %*% t(rep(1, K)))
      sigma <- 0
      for (i in 1:K) {
        for (j in 1:nrow(samp)) {
          sigma <- sigma + samp.wts[j] * dm.wts[j, i] *
                    sum((samp[j, ] - center[i, ])^2)
        }
      }
      sigma <- sqrt(sigma / p)
    }

    # --- Resampling (centre adaptation) ---
    center <- samp[resample(samp, K, samp.wts, method = resampling), ]
    center.all <- rbind(center.all, center)

    # --- Visualization ---
    if (visualization) {
      graphics::contour(x1, x2, z, drawlabels = FALSE, nlevels = 15,
                        main = sprintf("t = %d", t))
      graphics::points(samp, pch = 18, cex = 0.5, col = "green")
      graphics::points(center, pch = 16, cex = 1, col = "red")
    }
  }

  # --- Output ---
  if (output == "estimator") {
    z.std <- mean(exp(samp.all.logwts))
    samp.std.wts <- exp(samp.all.logwts - logaddexp(samp.all.logwts))
    m.std <- c(t(samp.all) %*% samp.std.wts)
    alpha <- ess / sum(ess)
    samp.wts.logwts <- samp.all.logwts + rep(log(alpha), each = K * J) + log(steps)
    z.wts <- mean(exp(samp.wts.logwts))
    samp.wts.wts <- exp(samp.wts.logwts - logaddexp(samp.wts.logwts))
    m.wts <- c(t(samp.all) %*% samp.wts.wts)
    m.las <- colMeans(center)
    return(list(m.std = m.std, m.wts = m.wts, m.las = m.las,
                z.std = z.std, z.wts = z.wts))
  } else {
    return(list(samp.all = samp.all,
                samp.all.logwts = samp.all.logwts,
                center.all = center.all,
                ess = ess))
  }
}
