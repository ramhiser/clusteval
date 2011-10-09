#' Ramey's Alternative Method
#'
#' TODO
#'
#' @export
#' @param x data matrix with N observations (rows) and p features (columns).
#' @param num_clusters TODO
#' @param K TODO
#' @param similarity_method TODO
#' @param noise TODO
#' @param B TODO
#' @param use_multicore TODO
#' @param ncpus TODO
#' @param ... TODO
#' @return list of scores by omitted cluster
clustomit_ramey <- function(x, K, cluster_method, similarity_method, B = 100, with_replacement = TRUE, use_multicore = FALSE, ncpus = 1, ...) {
  K <- as.integer(K)
  cluster_method <- as.character(cluster_method)
  similarity_method <- as.character(similarity_method)

  parallel_type <- ifelse(use_multicore, "multicore", "no")
  
  clusters <- cluster_wrapper(x, num_clusters = K, method = cluster_method, ...)
  out <- boot(x, ramey_boot, R = B, cluster_method = cluster_method,
    similarity_method = similarity_method, K = K, with_replacement = with_replacement,
    parallel = parallel_type, ncpus = ncpus, clusters = clusters
  )

	obj <- list(
		observed_value = out$t0,
		observed_mean = mean(out$t0),
		observed_min = min(out$t0),
		observed_maxn = max(out$t0),
		scores = out$t,
		mean = rowMeans(out$t),
		min = apply(out$t, 1, min),
		max = apply(out$t, 1, max),
		K = K,
		cluster_method = cluster_method,
		similarity_method = similarity_method,
		with_replacement = with_replacement,
		approach = "ramey"
	)
	class(obj) <- "clustomit"
	obj
}


#' The worker function for Ramey's Alternative Method
#'
#' This function is repeatedly called by boot() from the 'boot' package
#' to generate the approximate sampling distribution of Ramey's alternative method.
#'
#' TODO
#'
#' @param x data matrix with N observations (rows) and p features (columns).
#' @param idx TODO
#' @param K TODO
#' @param cluster_method TODO
#' @param similarity_method TODO
#' @param with_replacement TODO
#' @param ... TODO
#' @return list with results (TODO: Add more detail)
ramey_boot <- function(x, idx, K, cluster_method, similarity_method, with_replacement = TRUE, clusters, ...) {
  if(!with_replacement) {
    idx <- unique(idx)
  }
  idx <- sort(idx)
  cluster_levels <- unique(clusters)
  omit_similarities <- sapply(cluster_levels, function(k) {
    kept <- which(k != clusters)
    kept <- idx[idx %in% kept]
    x <- x[kept,]
    if(is.vector(x)) {
      x <- t(x)
    }
    
    if(nrow(x) == 0) {
      return(NA)
    }
    clusters_omit <- try(cluster_wrapper(x, num_clusters = K - 1, method = cluster_method, ...))
    if(inherits(clusters_omit, "try-error")) {
      warning("Returning NA")
      return(NA)
    }
    cluster_similarity(clusters[kept], clusters_omit, method = similarity_method)
  })
  # Unlike Landon's method, the cluster labelings are not arbitrary for each bootstrap replication
  # because the cluster labelings are set with the initial clustering of the data. Hence, we
  # are able to discuss the stability of individual clusters. Furthermore, similar to Landon's method
  # we can summarize each bootstrap replication (e.g. min, max, or mean) to better understand
  # the aggregate cluster stability of the current bootstrapped iteration.
  omit_similarities
}