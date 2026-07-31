\name{entropy}
\alias{entropy}
\title{Entropy Criterion}
\description{
Entropy (determinant of correlation matrix):
\deqn{\mathrm{Ent}(D)=|\boldsymbol{R}|}
}
\usage{
entropy(D, theta)
}
\arguments{
  \item{D}{Design points (numeric vector).}
  \item{theta}{Correlation length parameter.}
}
\value{Scalar determinant.}
