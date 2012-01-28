#' Cluster Data with a Selected Clustering Algorithm
#'
#' This function is a wrapper for a wide variety of clustering algorithms. For a
#' specified clustering algorithm, we apply the algorithm to a given data set
#' for a specified number of clusters
#'
#' The 'random' clustering algorithm randomly chooses one of K labels to apply to
#' each observation. By default, the probabilities for each cluster are equal,
#' but these can be customized.
#'
#' Custom algorithms are not yet allowed, but this feature is planned for a
#' future release.
#'
#' @param x a matrix containing the data to cluster. The rows are the sample
#' observations, and the columns are the features.
#' @param K the number of clusters
#' @param method a string specifying the clustering algorithm to apply. Possible values are: 'diana', 'kmeans', 'mclust', 'pam', and 'random'. See the Details section below.
#' @param ... optional arguments passed to the clustering algorithm
#' @return an object of class 'cluster_results'
cluster <- function(x, K, method = "kmeans", ...) {
  
}

cluster_lookup <- data.frame(name = "kmeans", package = "stats",
                             method = "kmeans", stringsAsFactors = FALSE) 
cluster_lookup <- rbind(cluster_lookup, c("diana", "cluster", "diana"))
cluster_lookup <- rbind(cluster_lookup, c("random", NULL, "random"))
cluster_lookup <- rbind(cluster_lookup, c("pam", "cluster", "pam"))
cluster_lookup <- rbind(cluster_lookup, c("model", "mclust", "Mclust"))
