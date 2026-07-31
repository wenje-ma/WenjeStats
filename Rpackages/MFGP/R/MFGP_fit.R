#' Fit a Multi-Fidelity Gaussian Process Model
#'
#' Implements the model of Tuo et al. (2014):
#' \deqn{y(\mathbf{x}, t) = \varphi(\mathbf{x}) + \delta(\mathbf{x}, t)}
#' where \eqn{\varphi(\mathbf{x})} is a stationary GP with anisotropic
#' Gaussian kernel, and \eqn{\delta(\mathbf{x}, t)} is a non-stationary GP
#' with the Gu et al. (2018) kernel on the fidelity parameter \eqn{t}.
#'
#' The covariance structure is:
#' \deqn{\mathrm{Cov}(y_1, y_2) = \sigma^2\big[R_0(\mathbf{x}_1,\mathbf{x}_2)
#'   + \lambda R_1(\mathbf{x}_1,\mathbf{x}_2) K(t_1, t_2)\big]}
#' where \eqn{R_0, R_1} are anisotropic Gaussian kernels and
#' \eqn{K(t_1, t_2) = e^{-(t_1-t_2)^2/\gamma^2} - e^{-t_1^2/\gamma^2}
#'   - e^{-t_2^2/\gamma^2} + 1}.
#'
#' Parameters are estimated via profile maximum likelihood using \code{\link{nlminb}}.
#'
#' @param D Design matrix. The first \code{p} columns are the input variables
#'   \eqn{\mathbf{x}}, and the last column is the fidelity parameter \eqn{t}.
#' @param y Response vector of length \code{nrow(D)}.
#' @param ini Initial values for the \code{2p+2} parameters: \eqn{\theta_0} (p),
#'   \eqn{\theta_1} (p), \eqn{\gamma} (1), \eqn{\lambda} (1).
#'   Default: \code{rep(1, 2*p+2)}.
#' @param nug Nugget term added to the diagonal of the correlation matrix for
#'   numerical stability. Default: \code{1e-6}.
#'
#' @return A list of class \code{"MFGP"} with components:
#'   \item{mu}{Estimated mean of the GP.}
#'   \item{nu2}{Estimated variance \eqn{\sigma^2}.}
#'   \item{theta}{Fitted parameters \eqn{(\theta_0, \theta_1, \gamma, \lambda)}.}
#'   \item{L}{Cholesky factor of the fitted correlation matrix.}
#'   \item{a}{Forward-solve of ones through \eqn{L^\top}.}
#'   \item{b}{Forward-solve of \eqn{y} through \eqn{L^\top}.}
#'   \item{D, dx2, t_vec, p}{Internal data for prediction.}
#'
#' @references
#' Tuo, R., Wu, C. F. J., & Yu, D. (2014). "Surrogate modeling of computer
#'   experiments with different mesh densities." \emph{Technometrics}, 56(3), 372-383.
#'
#' Gu, M., Wang, L., & Berger, J. O. (2018). "Nonstationary Gaussian process
#'   with application to computer emulation." \emph{Statistica Sinica}, 28(4), 2627-2647.
#'
#' @examples
#' # Simple 1D example with 2 input variables + fidelity
#' set.seed(1)
#' n <- 30
#' p <- 1
#' D <- cbind(runif(n), runif(n, 0.5, 1))
#' y <- sin(2 * pi * D[,1]) + 0.2 * (1 - D[,2]) * rnorm(n)
#' fit <- MFGP_fit(D, y)
#' pred <- MFGP_predict(fit, cbind(seq(0, 1, length.out = 20), 0))
#' plot(pred$mean, type = "l")
#'
#' @export
MFGP_fit <- function(D, y, ini = NULL, nug = 1e-6) {
  n <- length(y)
  if (!is.matrix(D)) D <- as.matrix(D)
  p <- ncol(D) - 1
  if (p < 1) stop("D must have at least 2 columns (x and t).")
  if (is.null(ini)) ini <- rep(1, 2 * p + 2)

  # Pre-compute per-dimension squared distance matrices
  dx2 <- lapply(1:p, function(k) outer(D[, k], D[, k], "-")^2)
  t_vec <- D[, p + 1]
  T1 <- matrix(t_vec, n, n)
  T2 <- matrix(t_vec, n, n, byrow = TRUE)
  Dt <- (T1 - T2)^2
  T1sq <- T1^2
  T2sq <- T2^2
  one <- rep(1, n)

  # Profile negative log-likelihood
  nLL <- function(theta) {
    s0 <- s1 <- 0
    for (k in 1:p) {
      s0 <- s0 + dx2[[k]] / theta[k]^2
      s1 <- s1 + dx2[[k]] / theta[p + k]^2
    }
    R0 <- exp(-s0)
    R1 <- exp(-s1)
    K <- exp(-Dt / theta[2 * p + 1]^2) -
         exp(-T1sq / theta[2 * p + 1]^2) -
         exp(-T2sq / theta[2 * p + 1]^2) + 1
    R <- R0 + theta[2 * p + 2] * R1 * K + diag(nug, n)
    L <- tryCatch(chol(R), error = function(e) return(NULL))
    if (is.null(L)) return(1e15)
    a <- forwardsolve(t(L), one)
    b <- forwardsolve(t(L), y)
    mu <- sum(a * b) / sum(a^2)
    nu2 <- max(1 / n * sum((b - mu * a)^2), 1e-15)
    log(nu2) + 2 * mean(log(diag(L)))
  }

  fit <- nlminb(ini, nLL,
                lower = ini / 10,
                upper = ini * 10,
                control = list(iter.max = 100))
  theta <- fit$par

  # Final MLE computation with fitted theta
  s0 <- s1 <- 0
  for (k in 1:p) {
    s0 <- s0 + dx2[[k]] / theta[k]^2
    s1 <- s1 + dx2[[k]] / theta[p + k]^2
  }
  R0 <- exp(-s0)
  R1 <- exp(-s1)
  K <- exp(-Dt / theta[2 * p + 1]^2) -
       exp(-T1sq / theta[2 * p + 1]^2) -
       exp(-T2sq / theta[2 * p + 1]^2) + 1
  R <- R0 + theta[2 * p + 2] * R1 * K + diag(nug, n)
  L <- chol(R)
  a <- forwardsolve(t(L), one)
  b <- forwardsolve(t(L), y)
  mu <- sum(a * b) / sum(a^2)
  nu2 <- max(1 / n * sum((b - mu * a)^2), 1e-15)

  structure(
    list(mu    = mu,
         nu2   = nu2,
         theta = theta,
         L     = L,
         a     = a,
         b     = b,
         D     = D,
         dx2   = dx2,
         t_vec = t_vec,
         p     = p),
    class = "MFGP"
  )
}
