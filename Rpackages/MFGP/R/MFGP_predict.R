#' Predict from a Fitted Multi-Fidelity GP Model
#'
#' Produces kriging predictions (mean and standard deviation) at new input
#' locations using a fitted model from \code{\link{MFGP_fit}}.
#'
#' @param obj A fitted model object of class \code{"MFGP"} returned by
#'   \code{\link{MFGP_fit}}.
#' @param new Matrix of new points. Must have \code{p + 1} columns, where
#'   \code{p} is the number of input variables. The first \code{p} columns
#'   are \eqn{\mathbf{x}} and the last column is the fidelity parameter
#'   \eqn{t}. To predict the true response (\eqn{t = 0}), set the last
#'   column to 0.
#'
#' @return A list with components:
#'   \item{mean}{Vector of predicted mean at each new point.}
#'   \item{sd}{Vector of kriging standard deviations.}
#'
#' @examples
#' set.seed(1)
#' n <- 30; p <- 1
#' D <- cbind(runif(n), runif(n, 0.5, 1))
#' y <- sin(2 * pi * D[,1]) + 0.2 * (1 - D[,2]) * rnorm(n)
#' fit <- MFGP_fit(D, y)
#' new <- cbind(seq(0, 1, length.out = 20), 0)
#' pred <- MFGP_predict(fit, new)
#' plot(new[,1], pred$mean, type = "l")
#' lines(new[,1], pred$mean + 2 * pred$sd, lty = 2, col = "gray")
#' lines(new[,1], pred$mean - 2 * pred$sd, lty = 2, col = "gray")
#'
#' @export
MFGP_predict <- function(obj, new) {
  if (!inherits(obj, "MFGP"))
    stop("obj must be of class 'MFGP' from MFGP_fit().")

  mu    <- obj$mu
  nu2   <- obj$nu2
  theta <- obj$theta
  L     <- obj$L
  a     <- obj$a
  b     <- obj$b
  D     <- obj$D
  dx2   <- obj$dx2
  t_vec <- obj$t_vec
  p     <- obj$p
  n     <- nrow(D)

  if (!is.matrix(new)) new <- as.matrix(new)
  if (ncol(new) != p + 1)
    stop(sprintf("new must have %d columns (p x-variables + t).", p + 1))

  theta0 <- theta[1:p]
  theta1 <- theta[(p + 1):(2 * p)]
  gam    <- theta[2 * p + 1]
  lam    <- theta[2 * p + 2]

  # Prediction at a single point
  ok <- function(x) {
    s0 <- s1 <- 0
    for (k in 1:p) {
      dk <- (D[, k] - x[k])^2
      s0 <- s0 + dk / theta0[k]^2
      s1 <- s1 + dk / theta1[k]^2
    }
    kt <- exp(-((t_vec - x[p + 1]) / gam)^2) -
          exp(-(t_vec / gam)^2) -
          exp(-(x[p + 1] / gam)^2) + 1
    r <- exp(-s0) + lam * exp(-s1) * kt
    d <- forwardsolve(t(L), r)
    pred_mu <- mu + sum(d * (b - mu * a))
    pred_sd <- sqrt(nu2 * max(1 - sum(d^2) + (1 - sum(d * a))^2 / sum(a^2), 0))
    c(pred_mu, pred_sd)
  }

  pred <- t(apply(new, 1, ok))
  list(mean = pred[, 1], sd = pred[, 2])
}
