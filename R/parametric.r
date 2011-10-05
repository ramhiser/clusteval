#' A Parametric Bootstrap Approach to Evaluate Clustering Stability via Cluster Omission.
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
clustomit_parboot <- function(x, num_clusters, cluster_method, similarity_method, noise = 0.1, B = 100, ...) {
	# Find the cluster for each observation in the provided data set, x.
	clust_orig <- cluster_wrapper(x, num_clusters, method = cluster_method, ...)

	# Iterates through each cluster and omit all of its observations for sensitivity analysis.
	clust_omit <- foreach(cl = levels(factor(clust_orig))) %do% {
		kept <- which(clust_orig != cl)
		x_kept <- x[kept,]
		clust_kept <- clust_orig[kept]

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

#' Perturbs a multivariate data set stored in a matrix.
#'
#' On the blog "X-Informatics" there is a post called "Add noise to data"
#' where a method is described to add noise to a data set. The method
#' finds the smallest dimension of the bounding box enclosing the data.
#' Take the square root of this, multiply by 'noise' (which is defaulted to 0.02), and
#' use this as the standard deviation when generating a normal random variate.
#'
#' The idea originates in the paper:
#' Balasubramanian, M., Schwartz, E.L., Tenenbaum, J.B., de Silva, V. & Langford, J.C.
#' "The Isomap Algorithm and Topological Stability". Science 295, 7a (2002).
#'
#' The blog post is located:
#' http://c13s.wordpress.com/2010/04/13/add-noise-to-data/
#'
#' A note about the minimum bounding box idea is located:
#' http://en.wikipedia.org/wiki/Minimum_bounding_box
#'
#' @export
#' @param x TODO
#' @param noise TODO
#' @return TODO
perturb_data <- function(x, noise = 0.1) {
	n <- nrow(x)
	p <- ncol(x)
	smallest_bounding_box_dim <- min(apply(x, 2, function(col) max(col) - min(col)))
	x + replicate(p, rnorm(n, sd = noise * sqrt(smallest_bounding_box_dim)))
}