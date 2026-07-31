#' Hilbert Curve Ordering
#'
#' Orders points by their distance along a Hilbert space-filling curve.
#' Useful for dimensionality reduction and stratification in resampling.
#'
#' @param x numeric matrix of points (rows = samples, columns = dimensions).
#' @param order integer; number of iterations of the Hilbert curve (default 8).
#'
#' @return integer vector of indices that sort \code{x} along the Hilbert curve.
#'
#' @importFrom gmp as.bigz as.bigq
#' @export
#'
#' @examples
#' x <- matrix(runif(20), ncol = 2)
#' hilbert_order(x)
hilbert_order <- function(x, order = 8){
  dim <- ncol(x)
  for (i in 1:dim){
    xi.min <- min(x[,i])
    xi.max <- max(x[,i])
    xi.range <- xi.max - xi.min
    xi.min <- xi.min - 1/(2^(order-1)) * xi.range
    xi.max <- xi.max + 1/(2^(order-1)) * xi.range
    x[,i] <- as.integer((x[,i] - xi.min) / (xi.max - xi.min) * 2^(order))
  }
  x.hdist <- apply(x, 1, hilbertcurve_distance_from_coordinates, dim = dim, order = order)
  x.order.idx <- order(x.hdist)
  return(x.order.idx)
}

# ---- Hilbert curve internals ----

# Convert decimal (bigz) to binary vector
DecToBin <- function(decimal, width = 64){
  if (!is.bigz(decimal)) stop("input must be a bigz object!")
  binary <- rep(0, width)
  n <- decimal
  i <- 1
  while (n > 0){
    r <- n %% 2
    binary[i] <- as.integer(r)
    n <- (n - r) / 2
    n <- as.bigz(n)
    i <- i + 1
  }
  binary <- rev(binary)
  return(binary)
}

# Convert binary vector to decimal (bigz)
BinToDec <- function(binary){
  if (any(!(binary %in% c(0,1)))) stop("input must be a binary object")
  p <- length(binary)
  binary <- rev(binary)
  decimal <- as.bigz(0)
  for (i in 1:p){
    decimal <- decimal + as.bigz(2^(i-1)) * binary[i]
  }
  return(decimal)
}

# Hilbert curve: distance → coordinates
hilbertcurve_coordinates_from_distance <- function(d, dim, order){
  if (d < 0 | d > 1) stop("distance must be normalised to 0 and 1!")
  if (dim < 1) stop("dimension must be positive!")
  if (order < 1) stop("order must be positive!")

  d <- d * as.bigq(as.bigz(2^(order * dim)))
  d <- as.bigz(d)
  d.bit <- DecToBin(d, order * dim)
  x <- rep(0, dim)
  for (i in 1:dim){
    x[i] <- as.integer(BinToDec(d.bit[seq(from = i, by = dim, length.out = order)]))
  }

  t <- bitwShiftR(x[dim], 1)
  for (i in dim:2) x[i] <- bitwXor(x[i], x[(i - 1)])
  x[1] <- bitwXor(x[1], t)

  Z <- bitwShiftL(2, (order - 1))
  Q <- 2
  while (Q != Z){
    P <- Q - 1
    for (i in dim:1){
      if (bitwAnd(x[i], Q)){
        x[1] <- bitwXor(x[1], P)
      } else {
        t <- bitwAnd(bitwXor(x[1], x[i]), P)
        x[1] <- bitwXor(x[1], t)
        x[i] <- bitwXor(x[i], t)
      }
    }
    Q <- bitwShiftL(Q, 1)
  }

  return(x)
}

# Hilbert curve: coordinates → distance
hilbertcurve_distance_from_coordinates <- function(x, dim, order){
  if (dim < 1) stop("dimension must be positive!")
  if (order < 1) stop("order must be positive!")
  if (length(x) != dim) stop(sprintf("x does not have dimension %d!", dim))
  if (any(x < 0)) stop("invalid coordinate: negative value!")
  if (any(x > (2^order - 1))) stop("invalid coordinate: exceeds 2^order - 1!")

  M <- bitwShiftL(1, order - 1)
  Q <- M
  while (Q > 1){
    P <- Q - 1
    for (i in 1:dim){
      if (bitwAnd(x[i], Q)){
        x[1] <- bitwXor(x[1], P)
      } else {
        t <- bitwAnd(bitwXor(x[1], x[i]), P)
        x[1] <- bitwXor(x[1], t)
        x[i] <- bitwXor(x[i], t)
      }
    }
    Q <- bitwShiftR(Q, 1)
  }

  for (i in 2:dim) x[i] <- bitwXor(x[i], x[(i - 1)])
  t <- 0
  Q <- M
  while (Q > 1){
    if (bitwAnd(x[dim], Q)) t <- bitwXor(t, Q - 1)
    Q <- bitwShiftR(Q, 1)
  }
  for (i in 1:dim) x[i] <- bitwXor(x[i], t)

  x.bit <- rep(0, order * dim)
  for (i in 1:dim){
    x.bit[seq(from = i, by = dim, length.out = order)] <- DecToBin(as.bigz(x[i]), order)
  }
  d <- BinToDec(x.bit)
  d <- as.numeric(d / as.bigz(2^(order * dim)))
  return(d)
}
