#' Consensus Clustering
#'
#' TODO
#' TODO: Add to the list of clustering algorithms
#'
#' @param x the data matrix
#' @param num_clusters the number of clusters
#' @param methods a vector of clustering methods to apply in consensus
#' @param ... (optional) additional arguments to pass to the clustering methods
#' @return TODO
#' @export
consensus <- function(x, num_clusters, method = c("hierarchical", "kmeans", "pam", "model", "diana"),
                      resample = c("subsample", "boot", "perturb"), num_resamples = 100,
                      subsample_size = 100, perturb_noise = 0.1, use_multicore = FALSE, ncpus = 1,
                      ...) {
  if (resample == "boot" || resample == "perturb") {
    stop("The bootstrap and perturbation methods are not implemented.")
  }
  # TODO: Use the boot function here with parallel. See clustomit for an example.
  # parallel_type <- ifelse(use_multicore, "multicore", "no")
  
  n <- nrow(x)
  
  out <- replicate(num_resamples, {
  subsample_idx <- sort(sample(x = seq.int(n), size = subsample_size))
  subsample_clust <- cluster_wrapper(x = x[subsample_idx, ],
                                     num_clusters = num_clusters,
                                     method = method)
  
  combs <- combn(subsample_idx, 2)
  list(
    comember_mat = spMatrix(nrow = n, ncol = n, i = combs[1, ], j = combs[2, ], x = comembership(subsample_clust)),
    indicator_mat = spMatrix(nrow = n, ncol = n, i = combs[1, ], j = combs[2, ], x = rep(1, ncol(combs)))
  )}, simplify = FALSE)

  comember_mat_sum <- Reduce("+", lapply(out, function(x) x$comember_mat))
  indicator_mat_sum <- Reduce("+", lapply(out, function(x) x$indicator_mat))
  consensus_mat <- as.matrix(comember_mat_sum / indicator_mat_sum)
  diag(consensus_mat) <- 0
  consensus_mat[which(is.nan(consensus_mat), arr.ind = TRUE)] <- 0

  # consensus_mat <- triu(comember_mat_sum / indicator_mat_sum, k = 1)
  list(comember_mat_sum = comember_mat_sum, indicator_mat_sum = indicator_mat_sum,
       consensus_mat = consensus_mat)
}

#' Heatmap of consensus matrix
#'
#' TODO
#' NOTE: Returns a ggplot2 object. Need to 'print' the object.
#'
#' @param mat the consensus matrix
#' @param num_clusters the number of clusters
#' @param methods a vector of clustering methods to apply in consensus
#' @param ... (optional) additional arguments to pass to the clustering methods
#' @return TODO
#' @export
plot_consensus <- function(mat) {
  require('reshape2')
  require('ggplot2')
  mat <- mat + t(mat)
  mat <- unname(mat)
  mat <- reshape2:::melt(mat)
  mat <- mat[complete.cases(mat), ]
  p <- ggplot(mat, aes(x=X1, y=X2, z = value, color = value))
  p <- p + geom_tile(aes(fill = value), colour = "white")
  p <- p + scale_fill_gradient(low = "white", high = "steelblue")
}

#' Randomly subsample a matrix or data frame.
#'
#' To randomly subsample, we are using the parametric option in the boot function
#' within the boot package. The 'parametric' option requires the specification of
#' a 'ran.gen' function to generate observations based on the original data and a list of
#' maximum likelihood estimators (MLE). We utilize this method, but instead of the
#' MLE, we instead pass the subsample_size.
#'
#' As noted in the boot documentation: Use of sim = "parametric" with a suitable ran.gen
#' allows the user to implement any types of nonparametric resampling which are not
#' supported directly.
#'
#' TODO: If 'x' and 'subsample_size' break, it may be necessary to use
#' 'data' and 'mle' instead.
#'
#' @param x the data matrix
#' @param subsample_size the number of observations (rows) to select at random from x.
#' @return a random subsample of the data matrix, x.
#' @export
#' @examples
#' TODO
boot_subsample <- function(x, subsample_size) {
  x[sample(x = seq_len(nrow(x)), size = subsample_size), ]
}

#' The worker function for consensus clustering
#' 
#' TODO
#'
#' @param x the data matrix
#' @param idx the row indices used with the bootstrap package, 'boot'.
#' @param num_clusters the number of clusters
#' @param methods a vector of clustering methods to apply in consensus
#' @param ... (optional) additional arguments to pass to the clustering methods
#' @return TODO
#' @export
consensus_boot <- function(x, idx, num_clusters, methods = c("hierarchical", "kmeans", "pam", "model", "diana"), ...) {
  x <- x[idx, ]
  cl_results <- llply(methods, cluster_wrapper, x = x, num_clusters = num_clusters)
}

#' A Wrapper Function for a List of Clustering Algorithms
#'
#' TODO
#' TODO: Add to the list of clustering algorithms
#' TODO: How to do everything in parallel?
#' TODO: For now, to keep things simple, I am using the already made clustering algorithm
#'       As such, renamed function to cluster_wrapper_mult. Eventually, use one and only one.
#'       This should be straightforward, but it would break too much for now.
#'
#' @param x the data matrix
#' @param K the number of clusters
#' @param methods a vector of clustering methods to apply in consensus
#' @param parallel If TRUE, applies each clustering algorithm in parallel. See Details.
#' @param ... (optional) additional arguments to pass to the clustering methods
cluster_wrapper_mult <- function(x, K, methods = c("kmeans", "pam", "hierarchical"), parallel = FALSE, ...) {
  warning("Not implemented yet.")
  NULL
}