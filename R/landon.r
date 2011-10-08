#' Landon's Alternative Method
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
clustomit_landon <- function(x, K, cluster_method, similarity_method, B = 100, with_replacement = TRUE, use_multicore = FALSE, ncpus = 1, ...) {
  K <- as.integer(K)
  cluster_method <- as.character(cluster_method)
  similarity_method <- as.character(similarity_method)
  
  parallel_type <- ifelse(use_multicore, "multicore", "no")
  
  out <- boot(x, landon_boot, R = B, cluster_method = cluster_method,
    similarity_method = similarity_method, K = K, with_replacement = with_replacement,
    parallel = parallel_type, ncpus = ncpus
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
		approach = "landon"
	)
	class(obj) <- "clustomit"
	obj
}


#' The worker function for Landon's Alternative Method
#'
#' This function is repeatedly called by boot() from the 'boot' package
#' to generate the approximate sampling distribution of Landon's alternative method.
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
landon_boot <- function(x, idx, K, cluster_method, similarity_method, with_replacement = TRUE, ...) {
  if(!with_replacement) {
    idx <- unique(idx)
  }
  x <- x[idx,]
  clusters <- cluster_wrapper(x, num_clusters = K, method = cluster_method, ...)
  cluster_levels <- unique(clusters)
  omit_similarities <- sapply(cluster_levels, function(k) {
    kept <- which(k != clusters)
    clusters_omit <- cluster_wrapper(x[kept,], num_clusters = K - 1, method = cluster_method, ...)
    cluster_similarity(clusters[kept], clusters_omit, method = similarity_method)
  })
  # With Landon's idea, it is important to understand that the ordering of 'omit_similarities'
  # will change each time because each call to 'landon_boot' yields different cluster labelings.
  # In other words, the cluster labelings are arbitrary for each bootstrap replication.
  # So, we can summarize each call to 'landon_boot' (e.g. min, max, or mean) to better understand
  # the cluster stability of the current bootstrapped iteration, but we cannot focus on the
  # stability of individual clusters without an ill-advised, computationally expensive ad-hoc approach,
  # which are often seen in the literature.
  omit_similarities
}