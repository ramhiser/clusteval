\name{boot_stratified_omit}
\alias{boot_stratified_omit}
\title{Creates a list of indices for a stratified nonparametric bootstrap.}
\usage{
  boot_stratified_omit(y, num_reps = 50)
}
\arguments{
  \item{y}{a vector that denotes the grouping of each
  observation. It must be coercible with \code{as.factor}.}

  \item{num_reps}{the number of bootstrap replications to
  use for each group}
}
\value{
  named list containing indices for each bootstrap
  replication
}
\description{
  This function creates a list of indices for a stratified
  nonparametric bootstrap. Corresponding to our Cluster
  Omission Stability statistic implemented in
  \code{\link{clustomit}}, we omit each group in turn and
  perform a stratified bootstrap without the group. We
  denote the number of groups as \code{num_clusters}, which
  is equal to \code{nlevels(factor(y))}. Specifically,
  suppose that we omit the \eqn{k}th group. That is, we
  ignore all of the observations corresponding to group
  \eqn{k}. Then, we sample with replacement from each of
  the remaining groups (i.e., every group except for group
  \eqn{k}), yielding a set of bootstrap indices.
}
\details{
  The returned list contains \eqn{K \times num_reps}
  elements.
}
\examples{
set.seed(42)
# We use 4 clusters, each with up to 10 observations. The sample sizes are
# randomly chosen.
num_clusters <- 4
sample_sizes <- sample(10, num_clusters, replace = TRUE)

# Create the cluster labels, y.
y <- unlist(sapply(seq_len(num_clusters), function(k) {
 rep(k, sample_sizes[k])
}))

# Use 20 reps per group.
boot_stratified_omit(y, num_reps = 20)

# Use the default number of reps per group.
boot_stratified_omit(y)
}

