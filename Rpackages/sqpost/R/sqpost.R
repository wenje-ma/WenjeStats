#' Gaussian Correlation Matrix
#'
#' Build the \eqn{n \times n} Gaussian correlation matrix for a set of
#' \eqn{p}-dimensional points, using dimension-specific length scales.
#'
#' @param x A numeric matrix of size \eqn{n \times p}.
#' @param s A numeric vector of length \eqn{p} giving the length scale
#'   for each dimension. Each element must be positive.
#'
#' @return An \eqn{n \times n} symmetric matrix with entries
#'   \deqn{G_{ij} = \exp\left\{-\frac12 \sum_{k=1}^p
#'          \frac{(x_{ik} - x_{jk})^2}{s_k^2}\right\}}
#'
#' @keywords internal
#' @examples
#' x <- matrix(runif(10), ncol = 2)
#' gauss_corr(x, s = c(0.5, 0.5))
gauss_corr <- function(x, s) {
  p <- ncol(x)
  # scaled distance matrix
  d <- as.matrix(dist(x %*% diag(1 / s, p, p), diag = TRUE, upper = TRUE))
  exp(-0.5 * d^2)
}


#' Estimate Length Scales via Cross-Validation
#'
#' Estimate the length scale parameters \eqn{s_1,\dots,s_p} for the
#' square-root posterior approximation by minimizing the log of the
#' weighted mean squared cross-validation error (wMSCV).
#'
#' The wMSCV for a candidate parameter vector \eqn{\phi = 1/s^2} is:
#' \deqn{\mathrm{wMSCV}(\phi) =
#'   \frac1n\sum_{i=1}^n \frac{[\boldsymbol{G}_\phi^{-1}
#'          \sqrt{\boldsymbol{h}}]_i^2}{[\boldsymbol{G}_\phi^{-1}]_{ii}}}
#' where \eqn{\boldsymbol{G}_\phi} is the Gaussian correlation matrix
#' computed with precision \eqn{\phi}.
#'
#' @param nu A numeric \eqn{n \times p} matrix of design points.
#' @param hval A numeric vector of length \eqn{n} giving the unnormalized
#'   posterior density (or any non-negative function) evaluated at each
#'   row of \code{nu}.
#' @param lambda A small positive regularization constant added to the
#'   diagonal of the correlation matrix for numerical stability.
#'   Default is \code{1e-4}.
#' @param lower,upper Optional numeric vectors of length \eqn{p} giving
#'   bounds for the initial guess. If \code{NULL}, they default to
#'   \code{ini.s/100} and \code{100*ini.s} respectively.
#' @param ... Additional arguments passed to \code{\link[stats]{optim}}.
#'
#' @return A numeric vector of length \eqn{p} containing the estimated
#'   length scale \eqn{\hat{s}_1,\dots,\hat{s}_p}.
#'
#' @export
#'
#' @references
#' Joseph, V. R. (2013). "Bayesian computation using design of experiments
#' and interpolation." \emph{Technometrics}, 55(4), 407-418.
#'
#' @examples
#' \dontrun{
#' set.seed(8)
#' library(MaxPro)
#' ini <- MaxPro(MaxProLHD(20, 2)$Design)$Design
#' library(mined)
#' nu <- mined(ini, logf, K_iter = 5)$cand
#' hev <- apply(nu, 1, h)
#' s <- estimate_scale(nu, hev)
#' }
estimate_scale <- function(nu, hval, lambda = 1e-4,
                           lower = NULL, upper = NULL, ...) {
  p <- ncol(nu)
  n <- nrow(nu)

  # Initial guess for phi = sigma^2 (variance)
  ini.phi <- (median(dist(nu)))^2 / 2

  # wMSCV objective: phi = sigma^2 (variance)
  wmscv <- function(phi) {
    G <- gauss_corr(nu, sqrt(phi))
    Gi <- solve(G + lambda * diag(n))
    mean(diag(Gi) * (c(Gi %*% sqrt(hval)) / diag(Gi))^2)
  }

  # Bounds for optimization (on phi = sigma^2 scale)
  if (is.null(lower)) lower <- rep(ini.phi / 100, p)
  if (is.null(upper)) upper <- rep(100 * ini.phi, p)

  opt <- optim(
    par = rep(ini.phi, p),
    fn  = function(phi) log(wmscv(phi)),
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    ...
  )

  # Return length scales s = sqrt(phi)
  sqrt(opt$par)
}


#' Fit Square-Root Posterior Approximation
#'
#' Fit the square-root basis expansion for an unnormalized posterior
#' density, following Joseph (2013). Given design points
#' \eqn{\boldsymbol{\nu}_1,\dots,\boldsymbol{\nu}_n} and evaluations
#' \eqn{h_i = h(\boldsymbol{\nu}_i)}, the approximation is
#' \deqn{\sqrt{h(\boldsymbol{\theta})} \approx
#'        \sum_{i=1}^n c_i \exp\!\left\{-\frac12\sum_{k=1}^p
#'          \frac{(\theta_k - \nu_{ik})^2}{s_k^2}\right\}}
#' where the coefficients \eqn{\boldsymbol{c}} are obtained by solving
#' \eqn{\boldsymbol{G}\boldsymbol{c} = \sqrt{\boldsymbol{h}}}.
#'
#' @param nu A numeric \eqn{n \times p} matrix of design points.
#' @param hval A numeric vector of length \eqn{n} with the unnormalized
#'   posterior evaluated at each design point.
#' @param s Optional numeric vector of length \eqn{p} giving the length
#'   scales. If \code{NULL} (default), they are estimated via
#'   \code{\link{estimate_scale}}.
#' @param lambda Regularization constant for the correlation matrix
#'   inverse. Default is \code{1e-4}.
#' @param ... Additional arguments passed to \code{\link{estimate_scale}}
#'   when \code{s} is \code{NULL}.
#'
#' @return An object of S3 class \code{"sqpost"} (invisibly), which
#'   is a list with components:
#'   \describe{
#'     \item{\code{nu}}{The \eqn{n \times p} design matrix.}
#'     \item{\code{hval}}{Evaluations of \eqn{h(\boldsymbol{\nu}_i)}.}
#'     \item{\code{s}}{Length scale vector \eqn{\hat{s}_1,\dots,\hat{s}_p}.}
#'     \item{\code{coef}}{Coefficient vector \eqn{\boldsymbol{c}}.}
#'     \item{\code{G}}{The \eqn{n \times n} Gaussian correlation matrix.}
#'     \item{\code{d}}{\eqn{n \times n} matrix with entries
#'           \eqn{d_{ij} = c_i c_j \exp\!\bigl(-\frac14\sum_k
#'             (\nu_{ik}-\nu_{jk})^2/s_k^2\bigr)}.}
#'     \item{\code{M}}{An \eqn{n \times n \times p} array where
#'           \eqn{M_{ijk} = (\nu_{ik} + \nu_{jk})/2}.}
#'   }
#'
#' @export
#'
#' @references
#' Joseph, V. R. (2013). "Bayesian computation using design of experiments
#' and interpolation." \emph{Technometrics}, 55(4), 407-418.
#'
#' @examples
#' \dontrun{
#' set.seed(8)
#' library(MaxPro)
#' ini <- MaxPro(MaxProLHD(20, 2)$Design)$Design
#' library(mined)
#' nu <- mined(ini, logf, K_iter = 5)$cand
#' hev <- apply(nu, 1, h)
#' fit <- sqrt_fit(nu, hev)
#' }
sqrt_fit <- function(nu, hval, s = NULL, lambda = 1e-4, ...) {
  p <- ncol(nu)
  n <- nrow(nu)

  # Step 1: Estimate or use provided length scales
  if (is.null(s)) {
    s <- estimate_scale(nu, hval, lambda = lambda, ...)
  }

  # Step 2: Build correlation matrix and solve for coefficients
  G <- gauss_corr(nu, s)
  G_reg <- G + lambda * diag(n)
  coef <- solve(G_reg) %*% sqrt(hval)

  # Step 3: Build pairwise matrices d and M
  # d_{ij} = c_i * c_j * exp(-0.25 * B^2_{ij})
  B <- as.matrix(dist(nu %*% diag(1 / s, p, p), diag = TRUE, upper = TRUE))
  d <- outer(c(coef), c(coef), "*") * exp(-0.25 * B^2)

  # M[,,k] = (nu_i_k + nu_j_k) / 2
  M <- array(0, dim = c(n, n, p))
  for (k in 1:p) {
    M[, , k] <- outer(nu[, k], nu[, k], "+") / 2
  }

  # Return S3 object
  structure(
    list(
      nu    = nu,
      hval  = hval,
      s     = s,
      coef  = coef,
      G     = G,
      d     = d,
      M     = M
    ),
    class = "sqpost"
  )
}


#' Marginal Posterior Density
#'
#' Compute the marginal posterior density for one dimension using the
#' fitted square-root approximation. The marginal density at a point
#' \eqn{\theta_k} is given by (Joseph, 2013, Eq. 4.53):
#' \deqn{\widehat{p}(\theta_k|\boldsymbol{y}) =
#'        \frac{\sum_{i=1}^n\sum_{j=1}^n d_{ij}\,
#'              \varphi\!\bigl(\theta_k; \frac{\nu_{ik}+\nu_{jk}}{2},
#'                              \frac{s_k^2}{2}\bigr)}
#'             {\sum_{i=1}^n\sum_{j=1}^n d_{ij}}}
#' where \eqn{\varphi} is the normal density.
#'
#' @param object An object of class \code{"sqpost"} returned by
#'   \code{\link{sqrt_fit}}.
#' @param grid A numeric vector of grid points at which to evaluate
#'   the marginal density.
#' @param dim Integer. The dimension (column index) for which to
#'   compute the marginal density. Default is 1.
#'
#' @return A numeric vector of the same length as \code{grid}, giving
#'   the marginal posterior density at each grid point.
#'
#' @export
#'
#' @references
#' Joseph, V. R. (2013). "Bayesian computation using design of experiments
#' and interpolation." \emph{Technometrics}, 55(4), 407-418.
#'
#' @examples
#' \dontrun{
#' fit <- sqrt_fit(nu, hev)
#' grid <- seq(0, 1, length.out = 100)
#' den1 <- marginal(fit, grid, dim = 1)
#' den2 <- marginal(fit, grid, dim = 2)
#' }
marginal <- function(object, grid, dim = 1) {
  stopifnot(inherits(object, "sqpost"))

  n  <- nrow(object$nu)
  p  <- ncol(object$nu)
  s  <- object$s[dim]
  d  <- object$d
  M_dim <- object$M[, , dim]

  # Density at each grid point via Eq. (4.53)
  den <- vapply(grid, function(th) {
    sum(d * dnorm(th, mean = M_dim, sd = s / sqrt(2)))
  }, numeric(1))

  # Normalize
  den / sum(d)
}
