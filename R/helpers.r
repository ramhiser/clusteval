#' A wrapper for the clustering wrappers
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param method TODO
#' @param ... TODO
#' @return TODO
cluster_wrapper <- function(x, num_clusters, method, ...) {
	cluster_methods <- c("kmeans", "singlelinkage", "pam", "diana", "mclust", "fanny")
	method <- match.arg(method, cluster_methods)
	method <- paste(method, "_wrapper", sep = "")
	get(method)(x = x, num_clusters = num_clusters, ...)
}

#' Wrapper for k-means, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param ... TODO
#' @return TODO
kmeans_wrapper <- function(x, num_clusters, ...) {
	kmeans(x = x, centers = num_clusters, ...)$cluster
}

#' Wrapper for single-linkage hierarhical, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param ... TODO
#' @return TODO
singlelinkage_wrapper <- function(x, num_clusters, ...) {
	cutree(hclust(dist(x), method = "single", ...), k = num_clusters)
}

#' Wrapper for PAM, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param ... TODO
#' @return TODO
pam_wrapper <- function(x, num_clusters, ...) {
	pam(x = x, k = num_clusters, ...)$cluster
}

#' Wrapper for DIANA, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param ... TODO
#' @return TODO
diana_wrapper <- function(x, num_clusters, ...) {
	cutree(as.hclust(diana(x = x, ...)), k = num_clusters)
}

#' Wrapper for model-based clustering, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param ... TODO
#' @return TODO
mclust_wrapper <- function(x, num_clusters, ...) {
	Mclust(data = x, G = num_clusters, ...)$classification
}


#' A wrapper for the similarity wrappers
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @param method TODO
#' @return TODO
cluster_similarity <- function(cl1, cl2, method) {
	similarity_measures <- c("rand", "jaccard")
	method <- match.arg(method, similarity_measures)
	method <- paste(method, "_wrapper", sep = "")
	
	# The clustering similarity measures prefer integer values for the cluster labels;
	# otherwise, a warning is thrown. The typecasting is done to avoid this warning.
	cl1 <- as.integer(cl1)
	cl2 <- as.integer(cl2)
	get(method)(cl1 = cl1, cl2 = cl2)
}

#' Wrapper for the Rand index for comparing two clusterings of the same data set.
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @return TODO
rand_wrapper <- function(cl1, cl2) {
	clv.Rand(std.ext(cl1, cl2))
}

#' Wrapper for the Jaccard index for comparing two clusterings of the same data set.
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @return TODO
jaccard_wrapper <- function(cl1, cl2) {
	clv.Jaccard(std.ext(cl1, cl2))
}