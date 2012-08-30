\name{comembership}
\alias{comembership}
\title{Calculates the comemberships of all pairs of a vector of clustering labels.}
\usage{
  comembership(labels)
}
\arguments{
  \item{labels}{a vector of \code{n} clustering labels}
}
\value{
  a vector of \code{choose(n, 2)} comembership bits
}
\description{
  For a set of clustering labels, this function computes
  the comembership of all pairs of observations. Basically,
  two observations are said to be comembers if they are
  clustered together.
}
\details{
  Tibshirani and Walther (2005) use the term
  'co-membership', which we shorten to 'comembership'. Some
  authors instead use the terms 'connectivity' or
  'co-occurrence'.

  We use the \code{Rcpp} package to improve the runtime
  speed of this function.
}
\examples{
# We generate K = 3 labels for each of n = 10 observations and compute the
# comembership for all 'n choose 2' pairs.
set.seed(42)
K <- 3
n <- 10
labels <- sample.int(K, n, replace = TRUE)
comembership_out <- comembership(labels)
comembership_out

# Notice that the number of comemberships is 'n choose 2'.
length(comembership_out) == choose(n, 2)
}
\references{
  Tibshirani, R. and Walther, G. (2005), Cluster Validation
  by Prediction Strength, _Journal of Computational and
  Graphical Statistics_, 14, 3, 511-528.
  \url{http://amstat.tandfonline.com/doi/abs/10.1198/106186005X59243}.
}

