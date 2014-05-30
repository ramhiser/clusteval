#' Creates a list of indices for a stratified nonparametric bootstrap.
#'
#' This function creates a list of indices for a nonparametric bootstrap.
#' Corresponding to our ClustOmit statistic implemented in
#' \code{\link{clustomit}}, we omit each cluster in turn and then sample from
#' the remaining clusters. We denote the number of groups as \code{K}, which is
#' equal to \code{nlevels(factor(y))}. Specifically, suppose that we omit the
#' \eqn{k}th group. That is, we ignore all of the observations corresponding to
#' group \eqn{k}. Then, we sample with replacement from each of the remaining
#' groups (i.e., every group except for group \eqn{k}), yielding a set of
#' bootstrap indices.
#'
#' The bootstrap resampling employed randomly samples from the remaining
#' observations after a cluster is omitted. By default, we ensure that one
#' observation is selected from each remaining cluster to avoid potential
#' situations where the resampled data set contains multiple replicates of a
#' single observation. Optionally, by setting the \code{stratified} argument to
#' \code{TRUE}, we employ a stratified sampling scheme, where instead we sample
#' with replacement from each cluster. In this case, the number of observations
#' sampled from a cluster is equal to the number of observations originally
#' assigned to that cluster (i.e., its cluster size).
#' The returned list contains \code{K * num_reps} elements.
#'
#' Both resampling schemes ensure that we avoid errors when clustering, similar
#' to this post on R Help:
#' \url{https://stat.ethz.ch/pipermail/r-help/2004-June/052357.html}.
#'
#' @export
#' @param y a vector that denotes the grouping of each observation. It must be
#' coercible with \code{\link{as.factor}}.
#' @param num_reps the number of bootstrap replications to use for each group
#' @return named list containing indices for each bootstrap replication
#' @param stratified Should the bootstrap replicates be stratified by cluster?
#' By default, no. See Details.
#' @examples
#' set.seed(42)
#' # We use 4 clusters, each with up to 10 observations. The sample sizes are
#' # randomly chosen.
#' K <- 4
#' sample_sizes <- sample(10, K, replace = TRUE)
#'
#' # Create the cluster labels, y.
#' y <- unlist(sapply(seq_len(K), function(k) {
#'  rep(k, sample_sizes[k])
#' }))
#'
#' # Use 20 reps per group.
#' boot_omit(y, num_reps = 20)
#'
#' # Use the default number of reps per group.
#' boot_omit(y)
boot_omit <- function(y, num_reps = 50, stratified = FALSE) {
  y <- as.factor(y)

  # Creates a list with each named element containing the indices of its cluster
  # members
  y_index <- split(seq_along(y), y)

  # We create a list, where each element contains the indices (rows) to apply
  # the clustering algorithm to examine its stability.
  boot_indices <- lapply(seq_len(nlevels(y)), function(omit_k) {
    y_omitted <- y_index[-omit_k]
    # If stratified is selected, then we apply 'sample_vector' to the indices
    # from each cluster.
    if (stratified) {
      replicate(n = num_reps, {
        unlist(lapply(y_omitted, sample_vector, replace = TRUE),
               use.names = FALSE)
      }, simplify = FALSE)
    } else {
      # By default, first we randomly sample a single observation from each
      # cluster. Then, we sample with replacement from the indices with the
      # current cluster omitted. The resulting number of labels is the length of
      # 'unlist(y_index[-omit_k])'
      replicate(n = num_reps, {
        first_draw <- unlist(lapply(y_omitted, sample_vector, size = 1),
               use.names = FALSE)
        y_omitted <- unlist(y_omitted, use.names = FALSE)
        second_draw <- sample(y_omitted, replace = TRUE,
                              size = length(y_omitted) - length(first_draw))
        c(first_draw, second_draw)
      }, simplify = FALSE)
    }
  })
  
  boot_indices <- unlist(boot_indices, recursive = FALSE)
  names(boot_indices) <- paste0("Rep", seq_along(boot_indices))
  boot_indices
}

#' Wrapper function to sample from a vector without any convenience features
#'
#' This function is a wrapper for the base \code{\link[base]{sample}} function
#' that circumvents the _undesired behavior_ from the _convenience feature_
#' discussed in the function's documentation.
#'
#' For example, if \code{x} contains a single element, say, 42, we wish
#' \code{sample(42)} to return only \code{42} rather than a random permutation
#' of \code{1:42}.
#' 
#' @param x vector
#' @param ... additional arguments passed to \code{\link[base]{sample}}
#' @return a vector containing a random permutation of \code{x}
sample_vector <- function(x, ...) {
  if (length(x) > 1) {
    sample(x, ...)
  } else {
    x
  }
}
