#' A Bootstrap Approach to Evaluate Clustering Stability via Cluster Omission.
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
#' @param cluster_method TODO
#' @param similarity_method TODO
#' @param noise TODO
#' @param B TODO
#' @param ... TODO
#' @return list of scores by omitted cluster
clustomit_boot <- function(x, num_clusters, cluster_method, similarity_method, B = 100, ...) {
	# Find the cluster for each observation in the provided data set, x.
	clust_orig <- cluster_wrapper(x, num_clusters, method = cluster_method, ...)
	
	lapply(clust_orig, function(k) {
	  kept <- which(clust_orig != cl)
		x_kept <- x[kept,]
		clust_kept <- clust_orig[kept]
	})

	# Iterates through each cluster and omit all of its observations for sensitivity analysis.
	clust_omit <- foreach(cl = levels(factor(clust_orig))) %do% {
		kept <- which(clust_orig != cl)
		x_kept <- x[kept,]
		clust_kept <- clust_orig[kept]
		
		boot()

		# Perturbs (jitters) the data B times and applies the user-specified clustering
		# method to the jittered data while omitting the current cluster of interest.
		clusters <- foreach(b=seq_len(B), .combine=rbind) %dopar% {
			x_perturbed <- perturb_data(x_kept, noise = noise)
			cluster_wrapper(x_perturbed, num_clusters - 1, method = cluster_method, ...)
		}
		
		# Compute the similarity of the original clustering with the clusterings from each perturbation.
		similarity <- apply(clusters, 1, function(clust) {
			cluster_similarity(cl1 = clust, cl2 = clust_kept, method = similarity_method)
		})
		list(clusters = clusters, similarity = as.vector(similarity), kept = kept)
	}
	
	list(
		clustering = factor(clust_orig),
		omit_results = clust_omit
	)
}

cluster_method <- "kmeans"
similarity_method <- "jaccard"
clustomit_worker <- function(x, idx, K, cluster_method, similarity_method, with_replacement = TRUE) {
  if(!with_replacement) {
    idx <- unique(idx)
  }
  x <- x[idx,]
  # TODO: Keep going
}