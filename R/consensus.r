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
consensus <- function(x, num_clusters, methods = c("hierarchical", "kmeans", "pam", "model", "diana"), ...) {
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