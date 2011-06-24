#' Evaluation of clustering stability by omitting each cluster at a time
#'
#' This is a bootstrapping approach to evaluating the stability of a clustering
#' algorithm through the clustering admissibility criteria in Fisher and 
#' Van Ness (1971). In particular, the focus here is on the cluster omission
#' admissibility criterion.
#'
#' We remove each cluster at a time. Then, we bootstrap from the given data
#' after removing the cluster, and use a measure of similarity between the
#' the cluster labels of the original data set and the cluster labels of
#' the bootstrapped sample.
#' @param x data matrix with N observations (rows) and p features (columns).
#' @param y vector of N labels corresponding to the rows of x 
#' @param cluster_f a reference to the clustering function being evaluated
#' @param ... additional arguments to pass to cluster_f
#' @return list of scores by omitted cluster
clust_omit <- function(x, y, cluster_f, ...) {
  clusters <- unique(y)

}
