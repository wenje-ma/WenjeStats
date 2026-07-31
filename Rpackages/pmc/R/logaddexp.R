#' Log-Sum-Exp (Numerically Stable)
#'
#' Computes \code{log(sum(exp(logv)))} in a numerically stable way,
#' avoiding underflow when \code{logv} contains large negative values.
#'
#' @param logv numeric vector of log-scale values.
#'
#' @return \code{log(sum(exp(logv)))}.
#'
#' @export
logaddexp <- function(logv){
  logv.max <- max(logv)
  logv.sum <- log(sum(exp(logv - logv.max))) + logv.max
  return(logv.sum)
}
