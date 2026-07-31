#' Energy Distance
#'
#' Computes the energy distance between two samples \code{X} and \code{Y},
#' optionally with weighted samples.
#'
#' The energy distance is defined as
#' \deqn{2E\|X-Y\| - E\|X-X'\|,}
#' where \eqn{X'} is an independent copy of \eqn{X}.
#'
#' @param X numeric matrix of samples from distribution P.
#' @param Y numeric matrix of samples from distribution Q.
#' @param X.wts optional numeric vector of weights for rows of \code{X}.
#' @param Y.wts optional numeric vector of weights for rows of \code{Y}.
#'
#' @return the energy distance (scalar).
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 2)
#' Y <- matrix(rnorm(100, mean = 1), ncol = 2)
#' energy_distance(X, Y)
energy_distance <- function(X, Y, X.wts = NULL, Y.wts = NULL) {
  n.X <- nrow(X)
  n.Y <- nrow(Y)
  if (is.null(X.wts)) X.wts <- rep(1 / n.X, n.X)
  if (is.null(Y.wts)) Y.wts <- rep(1 / n.Y, n.Y)
  e.dist <- 0
  for (i in 1:n.X) {
    for (j in 1:n.Y) {
      e.dist <- e.dist + 2 * X.wts[i] * Y.wts[j] * sqrt(sum((X[i, ] - Y[j, ])^2))
    }
  }
  for (i in 1:n.X) {
    for (j in 1:n.X) {
      e.dist <- e.dist - X.wts[i] * X.wts[j] * sqrt(sum((X[i, ] - X[j, ])^2))
    }
  }
  return(e.dist)
}
