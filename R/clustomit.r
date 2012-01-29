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
  weighted_mean = TRUE, B = 100, with_replacement = FALSE, use_multicore = FALSE,
  ncpus = 1, ...) {

  # TODO: Unit tests to make sure these arguments are handled correctly.
  K <- as.integer(K)
  cluster_method <- as.character(cluster_method)
  similarity_method <- as.character(similarity_method)

  parallel_type <- ifelse(use_multicore, "multicore", "no")

  # The cluster labels for the observed (original) data (i.e. x).  
  obs_clusters <- cluster_wrapper(x, num_clusters = K, method = cluster_method, ...)
  obs_num_clusters <- as.vector(table(obs_clusters))

  boot_out <- boot(x, clustomit_boot, R = B, cluster_method = cluster_method,
    similarity_method = similarity_method, K = K, with_replacement = with_replacement,
    parallel = parallel_type, ncpus = ncpus, obs_clusters = obs_clusters
  )

  if (weighted_mean) {
    obs_aggregate <- weighted.mean(boot_out$t0, w = obs_num_clusters)
    boot_aggregate <- apply(boot_out$t, 1, weighted.mean, w = obs_num_clusters)    
  } else {
    obs_aggregate <- mean(boot_out$t0)
    boot_aggregate <- rowMeans(boot_out$t)
  }

	obj <- list(
		obs_cluster_similarity = boot_out$t0,
		obs_aggregate = obs_aggregate,
		boot_cluster_similarity = boot_out$t,
    boot_aggregate = boot_aggregate,
    obs_clusters = obs_clusters,
		K = K,
		cluster_method = cluster_method,
		similarity_method = similarity_method,
		with_replacement = with_replacement
	)
	class(obj) <- "clustomit"
	obj
}


#' The worker function for Cluster Omission
#'
#' This function is repeatedly called by boot() from the 'boot' package
#' to generate the approximate sampling distribution of the similarity
#' index of the cluster omission method.
#'
#' TODO
#'
#' @param x data matrix with N observations (rows) and p features (columns).
#' @param idx TODO
#' @param K TODO
#' @param cluster_method TODO
#' @param similarity_method TODO
#' @param with_replacement TODO
#' @param obs_clusters vector cluster labels for each observation from the observed (original) data set
#' @param ... TODO
#' @return vector the similarity scores for each omitted cluster
clustomit_boot <- function(x, idx, K, cluster_method, similarity_method, with_replacement = FALSE, obs_clusters, ...) {
  if(!with_replacement) {
    idx <- unique(idx)
  }
  idx <- sort(idx)
  cluster_levels <- unique(obs_clusters)
  omit_similarities <- sapply(cluster_levels, function(k) {
    kept <- which(k != obs_clusters)
    kept <- idx[idx %in% kept]
    x <- x[kept,]
    if(is.vector(x)) {
      x <- t(x)
    }
    
    if(nrow(x) == 0) {
      return(NA)
    }
    clusters_omit <- try(cluster_wrapper(x, num_clusters = K - 1, method = cluster_method, ...),
                         silent = TRUE)

    if(inherits(clusters_omit, "try-error")) {
      warning("Error in clustomit_boot: ", attr(clusters_omit, "condition")$message, "\nNumber of Rows in x: ", nrow(x), "\nNumber of Clusters: ", K - 1)
      return(clusters_omit)
    }
    cluster_similarity(obs_clusters[kept], clusters_omit, method = similarity_method)
  })

  omit_similarities
}
