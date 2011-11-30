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
consensus <- function(x, num_clusters, methods = c("hierarchical", "kmeans", "pam", "model", "diana"),
                      resample = c("boot", "subsample", "perturb"), num_resamples = 100,
                      subsample_size = 100, perturb_noise = 0.1, use_multicore = FALSE, ncpus = 1,
                      ...) {
  if (resample == "boot" || resample == "perturb") {
    stop("The bootstrap and perturbation methods are not implemented.")
  }
  # TODO: Use the boot function here with parallel. See clustomit for an example.
  parallel_type <- ifelse(use_multicore, "multicore", "no")
  
  boot_out <- boot(x, # TODO: Statistic,
                   sim = 'parametric', R = num_resamples, methods = methods,
                   ran.gen = boot_subsample, mle = subsample_size,
                   parallel = parallel_type, ncpus = ncpus, ...)
  names(cl_results) <- methods
  cl_results
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
  x[sample(x = seq_len(n), size = subsample_size), ]
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