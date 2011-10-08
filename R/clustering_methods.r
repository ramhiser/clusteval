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
	cluster_methods <- c("hierarchical", "kmeans", "pam", "model", "diana")
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

#' Wrapper for hierarhical clustering, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' The default is single-linkage.
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param ... TODO
#' @return TODO
hierarchical_wrapper <- function(x, num_clusters, method = "single", ...) {
	cutree(hclust(dist(x), method = method, ...), k = num_clusters)
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
diana_wrapper <- function(x, num_clusters, hennig = FALSE, ...) {
	diana_out <- cutree(as.hclust(diana(x = x, ...)), k = num_clusters)
	if(hennig) {
    hennig_list <- list()
    hennig_list$result <- diana_out
    hennig_list$nc <- num_clusters
    hennig_list$partition <- diana_out
    hennig_list$clustermethod <- "diana"
    # A list whose elements are logical vectors of length n.
    # The ith element in the list is 1 if observation i is said to belong to the current class,
    # and 0 otherwise.
    hennig_list$clusterlist <- lapply(unique(diana_out), function(cl) {
      as.integer(cl == diana_out)
    })
    return(hennig_list)
	} else {
	 return(diana_out)
	}
}

#' Wrapper for model-based clustering, so that only the data set and number of clusters is specified.
#'
#' TODO
#'
#' To use Hennig's (2007) clusterboot function from the 'fpc' package with this wrapper for model-based
#' clustering, we have to return a custom list to 
#'
#' @export
#' @param x TODO
#' @param num_clusters TODO
#' @param hennig TODO
#' @param ... TODO
#' @return TODO
model_wrapper <- function(x, num_clusters, hennig = FALSE, ...) {
  mclust_out <- Mclust(data = x, G = num_clusters, ...)
  if(hennig) {
    hennig_list <- list()
    hennig_list$result <- mclust_out
    hennig_list$nc <- num_clusters
    hennig_list$partition <- mclust_out$classification
    hennig_list$clustermethod <- "mclust"
    # A list whose elements are logical vectors of length n.
    # The ith element in the list is 1 if observation i is said to belong to the current class,
    # and 0 otherwise.
    hennig_list$clusterlist <- lapply(unique(mclust_out$classification), function(cl) {
      as.integer(cl == mclust_out$classification)
    })
    return(hennig_list)
	} else {
	 return(mclust_out$classification)
	}
}

