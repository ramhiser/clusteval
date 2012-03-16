#' Cluster Stability Evaluation via Cluster Omission
#'
#' TODO: Add more thorough documentation
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
#'
#' The cluster method should take the data as 'x' and the 'num_clusters'.
#'
#' @export
#' @param x data matrix with N observations (rows) and p features (columns).
#' @param num_clusters TODO
#' @param K TODO
#' @param similarity_method TODO
#' @param weighted_mean logical Should the aggergate similarity score for each
#'  bootstrap replication be weighted by the number of observations in each of
#'  the observed clusters?
#' @param B TODO
#' @param use_multicore TODO
#' @param ncpus TODO
#' @param ... TODO
#' @return list of scores by omitted cluster
clustomit <- function(x, K, cluster_method, similarity_method = "jaccard",
                      weighted_mean = TRUE, B = 100, parallel = FALSE,
                      num_cores = getOption("mc.cores", 2), num_reps = 50, ...) {

  K <- as.integer(K)
  cluster_method <- as.character(cluster_method)
  similarity_method <- as.character(similarity_method)

  # If the user indicates that parallel should not be used, then set the number
  # of cores to 1. This forces 'mclapply' to use 'lapply.
  if (!parallel) {
    num_cores <- 1
  }

  # The cluster labels for the observed (original) data (i.e. x).  
  obs_clusters <- cluster_wrapper(x, num_clusters = K, method = cluster_method,
                                  ...)
  cluster_sizes <- as.vector(table(obs_clusters))

  # Determines the indices for the bootstrap reps.
  boot_indices <- boot_stratified_omit(obs_clusters, num_reps = num_reps)

  # For each set of bootstrap indices, cluster the resampled data and
  # compute the similarity with the corresponding original clusters.
  boot_similarity <- mclapply(boot_indices, function(idx) {
    clusters_omit <- cluster_wrapper(x[idx, ], num_clusters = K - 1,
                                     method = cluster_method, ...)

    cluster_similarity(obs_clusters[idx], clusters_omit,
                       method = similarity_method)
  }, mc.cores = num_cores)

  # Because 'mclapply' returns a list, we first simplify the list to an array and
  # then 'split' the similarity scores into a list by cluster.
  boot_similarity <- simplify2array(boot_similarity)
  boot_similarity <- split(boot_similarity, gl(K, num_reps))

  # Now, we compute the weighted average of the similarity scores for each
  # bootstrap replication. The weights correspond to the sample sizes of each
  # cluster.
  boot_similarity_matrix <- do.call(cbind, boot_similarity)
  boot_aggregate <- apply(boot_similarity_matrix, 1, weighted.mean,
                          w = cluster_sizes)

	obj <- list(
		boot_similarity = boot_similarity,
    boot_aggregate = boot_aggregate,
    obs_clusters = obs_clusters,
		K = K,
		cluster_method = cluster_method,
		similarity_method = similarity_method
	)
	class(obj) <- "clustomit"
	obj
}
