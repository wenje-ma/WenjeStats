#' Resample Weighted Samples
#'
#' Performs resampling of weighted samples using one of several methods.
#' This is a unified interface to multinomial, residual, systematic,
#' stratified, and support points resampling.
#'
#' @param x numeric matrix of samples (rows = samples, columns = dimensions).
#' @param n integer; number of resampled indices to return.
#' @param prob numeric vector of normalised weights (must sum to 1).
#' @param method character; one of \code{"multinomial"}, \code{"residual"},
#'        \code{"systematic"}, \code{"stratified"}, or \code{"sp"}.
#' @param tol tolerance for support points iterative refinement (default 1e-6).
#' @param iter.max maximum iterations for support points refinement (default 10).
#'
#' @return integer vector of indices (with replacement) into \code{x}.
#'
#' @export
#'
#' @examples
#' x <- matrix(rnorm(100), ncol = 2)
#' w <- rep(1/50, 50)
#' idx <- resample(x, 20, w, method = "multinomial")
#' idx <- resample(x, 20, w, method = "systematic")
#' idx <- resample(x, 20, w, method = "sp")
resample <- function(x, n, prob,
                     method = c("multinomial", "residual", "systematic",
                                "stratified", "sp"),
                     tol = 1e-6, iter.max = 10) {
  method <- match.arg(method)
  idx <- switch(method,
    multinomial = .resample_multinomial(x, n, prob),
    residual    = .resample_residual(x, n, prob),
    systematic  = .resample_systematic(x, n, prob),
    stratified  = .resample_stratified(x, n, prob),
    sp          = .resample_sp(x, n, prob, tol, iter.max)
  )
  return(idx)
}

# ---- Internal resampling implementations ----

.resample_multinomial <- function(x, n, prob) {
  N <- nrow(x)
  sample(1:N, n, replace = TRUE, prob = prob)
}

.resample_residual <- function(x, n, prob) {
  N <- nrow(x)
  ept <- n * prob
  cnt <- floor(ept)
  wts <- ept - cnt
  idx <- c()
  if (sum(wts) > 0) {
    wts <- wts / sum(wts)
    idx <- sample(1:N, n - sum(cnt), replace = TRUE, prob = wts)
  }
  for (i in 1:N) idx <- c(idx, rep(i, cnt[i]))
  return(idx)
}

.resample_systematic <- function(x, n, prob) {
  N <- nrow(x)
  x.order.idx <- hilbert_order(x)
  prob <- prob[x.order.idx]
  idx <- c()
  u <- runif(1) / n
  l <- 0
  j <- 0
  while (u < 1) {
    if (l > u) {
      u <- u + 1 / n
      idx <- c(idx, j)
    } else {
      j <- j + 1
      l <- l + prob[j]
    }
  }
  return(x.order.idx[idx])
}

.resample_stratified <- function(x, n, prob) {
  N <- nrow(x)
  x.order.idx <- hilbert_order(x)
  prob <- prob[x.order.idx]
  idx <- c()
  prob.cumsum <- cumsum(prob)
  endpoints <- seq(0, 1, length.out = n + 1)
  for (i in 1:n) {
    lb <- sum(prob.cumsum <= endpoints[i]) + 1
    ub <- min(sum(prob.cumsum <= endpoints[i + 1]) + 1, N)
    if (lb == ub) {
      idx <- c(idx, lb)
    } else {
      prob.loc <- prob[lb:ub]
      prob.loc[1] <- prob.cumsum[lb] - endpoints[i]
      prob.loc[length(prob.loc)] <- endpoints[i + 1] - prob.cumsum[ub - 1]
      prob.loc <- prob.loc / sum(prob.loc)
      idx <- c(idx, sample(lb:ub, 1, replace = TRUE, prob = prob.loc))
    }
  }
  return(x.order.idx[idx])
}

.resample_sp <- function(x, n, prob, tol, iter.max) {
  if (is.null(dim(x))) stop("x must be a matrix!")
  N <- nrow(x)
  x.dist <- as.matrix(dist(x))
  x.measure <- c(x.dist %*% prob)
  idx <- which.min(x.measure)
  if (prob[idx[1]] < 1 / N) idx[1] <- which.max(prob)
  for (i in 2:n) {
    measure <- x.measure - apply(matrix(x.dist[, idx], ncol = i - 1), 1, sum) / i
    idx <- c(idx, which.min(measure))
  }
  edist <- 2 * sum(x.measure[idx]) / n - sum(x.dist[idx, idx]) / n^2
  iter <- 0
  while (TRUE) {
    iter <- iter + 1
    for (i in 1:n) {
      measure <- x.measure - apply(x.dist[, idx[-i]], 1, sum) / n
      idx[i] <- which.min(measure)
    }
    edist.new <- 2 * sum(x.measure[idx]) / n - sum(x.dist[idx, idx]) / n^2
    if ((edist - edist.new) < tol || iter > iter.max) break
    edist <- edist.new
  }
  return(idx)
}
