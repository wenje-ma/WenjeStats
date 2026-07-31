#' Gaussian Correlation Function (1D)
#'
#' Compute the \eqn{n \times m} Gaussian correlation matrix between
#' two sets of 1D points.
#' \deqn{R(x_i, y_j) = \exp\left\{-\left(\frac{x_i - y_j}{\theta}\right)^2\right\}}
#'
#' @param x A numeric vector of length \eqn{n}.
#' @param y A numeric vector of length \eqn{m}. If \code{NULL} (default),
#'   uses \code{x} as \code{y}.
#' @param theta Positive numeric. Correlation length parameter.
#'
#' @return An \eqn{n \times m} matrix of correlations.
#' @export
#'
#' @examples
#' gauss_corr(0:5/5, theta = 0.2)
#' gauss_corr(0:5/5, 0:10/10, theta = 0.1)
gauss_corr <- function(x, y = NULL, theta) {
  if (is.null(y)) y <- x
  exp(-(outer(x, y, "-") / theta)^2)
}


#' GP RMSE (Root Mean Squared Error)
#'
#' Compute the standardized RMSE \eqn{s(\boldsymbol{x}) = \sqrt{1 - \boldsymbol{r}(\boldsymbol{x})'\boldsymbol{R}^{-1}\boldsymbol{r}(\boldsymbol{x})}}
#' at a set of test points, given a design \eqn{D} and correlation
#' parameter \eqn{\theta}. A small regularization is added to the
#' diagonal of \eqn{\boldsymbol{R}} for numerical stability.
#'
#' @param D Numeric vector of design points (length \eqn{n}).
#' @param theta Positive numeric. Correlation length parameter.
#' @param test Numeric vector of test points.
#' @param lambda Small regularization for \eqn{\boldsymbol{R}} inverse.
#'   Default \code{1e-6}.
#'
#' @return A numeric vector of RMSE values at \code{test}.
#' @export
#'
#' @examples
#' D <- (0:9) / 9
#' test <- seq(0, 1, length = 101)
#' rmse <- gp_rmse(D, theta = 0.1, test)
gp_rmse <- function(D, theta, test, lambda = 1e-6) {
  n <- length(D)
  R <- gauss_corr(D, theta = theta)
  R_inv <- solve(R + lambda * diag(n))
  r_test <- gauss_corr(test, D, theta = theta)
  Z <- r_test %*% R_inv
  sqrt(1 - rowSums(r_test * Z))
}


#' IMSE Criterion
#'
#' Compute the Integrated Mean Squared Error criterion for a design
#' \eqn{D} with correlation parameter \eqn{\theta}:
#' \deqn{\mathrm{IMSE}(D;\theta) = \int_0^1 (1 - \boldsymbol{r}(x)'\boldsymbol{R}^{-1}\boldsymbol{r}(x)) \,\mathrm{d}x}
#' The integral is evaluated via numerical integration over \eqn{[0,1]}.
#'
#' @param D Numeric vector of design points (length \eqn{n}).
#' @param theta Positive numeric. Correlation length parameter.
#' @param lambda Small regularization for \eqn{\boldsymbol{R}} inverse.
#'   Default \code{1e-6}.
#'
#' @return A scalar: the IMSE value.
#' @export
#'
#' @examples
#' D <- (0:9) / 9
#' imse(D, theta = 0.1)
imse <- function(D, theta, lambda = 1e-6) {
  n <- length(D)
  R <- gauss_corr(D, theta = theta)
  R_inv <- solve(R + lambda * diag(n))

  # Integrand: 1 - r(x)' R^{-1} r(x)
  fn <- function(x) {
    r_x <- gauss_corr(x, D, theta = theta)
    Z <- r_x %*% R_inv
    1 - sum(r_x * Z)
  }

  integrate(Vectorize(fn), 0, 1)$value
}


#' MMSE Criterion
#'
#' Compute the Maximum Mean Squared Error criterion for a design
#' \eqn{D} with correlation parameter \eqn{\theta}:
#' \deqn{\mathrm{MMSE}(D;\theta) = \max_{x \in [0,1]} (1 - \boldsymbol{r}(x)'\boldsymbol{R}^{-1}\boldsymbol{r}(x))}
#' The maximum is approximated by evaluating on a fine grid of test points.
#'
#' @param D Numeric vector of design points (length \eqn{n}).
#' @param theta Positive numeric. Correlation length parameter.
#' @param n_grid Integer. Number of grid points for approximating the
#'   maximum. Default \code{301}.
#' @param lambda Small regularization for \eqn{\boldsymbol{R}} inverse.
#'   Default \code{1e-6}.
#'
#' @return A scalar: the MMSE value.
#' @export
#'
#' @examples
#' D <- (0:9) / 9
#' mmse(D, theta = 0.1)
mmse <- function(D, theta, n_grid = 301, lambda = 1e-6) {
  test <- seq(0, 1, length = n_grid)
  max(gp_rmse(D, theta, test, lambda)^2)
}


#' Entropy Criterion
#'
#' Compute the entropy criterion for a design \eqn{D} with correlation
#' parameter \eqn{\theta}:
#' \deqn{\mathrm{Ent}(D) = |\boldsymbol{R}|}
#' where \eqn{\boldsymbol{R}_{ij} = \exp(-(D_i - D_j)^2/\theta^2)}.
#'
#' @param D Numeric vector of design points (length \eqn{n}).
#' @param theta Positive numeric. Correlation length parameter.
#'
#' @return A scalar: the determinant of the correlation matrix.
#' @export
#'
#' @examples
#' D <- (0:9) / 9
#' entropy(D, theta = 0.1)
entropy <- function(D, theta) {
  R <- gauss_corr(D, theta = theta)
  det(R)
}
