\name{boot_omit}
\alias{boot_omit}
\title{Creates a list of indices for a stratified nonparametric bootstrap.}
\usage{
  boot_omit(y, num_reps = 50, stratified = FALSE)
}
\arguments{
  \item{y}{a vector that denotes the grouping of each
  observation. It must be coercible with
  \code{\link{as.factor}}.}

  \item{num_reps}{the number of bootstrap replications to
  use for each group}
}
\value{
  named list containing indices for each bootstrap
  replication
}
\description{
  This function creates a list of indices for a
  nonparametric bootstrap. Corresponding to our ClustOmit
  statistic implemented in \code{\link{clustomit}}, we omit
  each cluster in turn and then sample from the remaining
  clusters. We denote the number of groups as \code{K},
  which is equal to \code{nlevels(factor(y))}.
  Specifically, suppose that we omit the \eqn{k}th group.
  That is, we ignore all of the observations corresponding
  to group \eqn{k}. Then, we sample with replacement from
  each of the remaining groups (i.e., every group except
  for group \eqn{k}), yielding a set of bootstrap indices.
}
\details{
  The bootstrap resampling employed randomly samples from
  the remaining observations after a cluster is omitted. By
  default, we ensure that one observation is selected from
  each remaining cluster to avoid potential situations
  where the resampled data set contains multiple replicates
  of a single observation. Optionally, by setting the
  \code{stratified} argument to \code{TRUE}, we employ a
  stratified sampling scheme, where instead we sample with
  replacement from each cluster. In this case, the number
  of observations sampled from a cluster is equal to the
  number of observations originally assigned to that
  cluster (i.e., its cluster size). The returned list
  contains \code{K * num_reps} elements.

  Both resampling schemes ensure that we avoid errors when
  clustering, similar to this post on R Help:
  \url{https://stat.ethz.ch/pipermail/r-help/2004-June/052357.html}.
}
\examples{
set.seed(42)
# We use 4 clusters, each with up to 10 observations. The sample sizes are
# randomly chosen.
K <- 4
sample_sizes <- sample(10, K, replace = TRUE)

# Create the cluster labels, y.
y <- unlist(sapply(seq_len(K), function(k) {
 rep(k, sample_sizes[k])
}))

# Use 20 reps per group.
boot_omit(y, num_reps = 20)

# Use the default number of reps per group.
boot_omit(y)
}

