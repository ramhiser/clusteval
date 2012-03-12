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

  # Several of the clustering algorithms throw an error if the number of
  # observations are equal to the number of clusters. In this case, each of the
  # clusters are singletons. So, we return a sequence of integers with the same
  # length as 'num_clusters'.
  if (nrow(x) == num_clusters) {
    return(seq_len(num_clusters))
  }
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
kmeans_wrapper <- function(x, num_clusters, num_starts = 10, ...) {
	kmeans(x = x, centers = num_clusters, nstart = num_starts, ...)$cluster
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


#' Randomly cluster a data set into K clusters.
#'
#' For each observation (row) in 'x', one of K labels is randomly generated.
#' By default, the probabilities of selecting each clustering label are equal,
#' but this can be altered by specifying 'prob', a vector of probabilities for
#' each cluster.
#'
#' Random clustering is often utilized as a baseline comparison clustering
#' against which other clustering algorithms are employed to identify structure
#' within the data. Typically, comparisons are made in terms of proposed
#' clustering assessment and evaluation methods as well as clustering similarity
#' measures. For the former, a specified clustering evaluation method is computed
#' for the considered clustering algorithms as well as random clustering. If the
#' clusters determined by a considered clustering algorithm do not differ
#' significantly from the random clustering, we might conclude that the found
#' clusters are no better than naively choosing clustering labels for each
#' observation at random. Likewise, a similarity measure can be computed to
#' compare the clusterings from each of a considered clustering algorithm and a
#' random clustering: if the clusterings are significantly similar, once again,
#' we might conclude the clusters found via the considered clustering algorithm
#' do not differ significantly from those found at random. In either case, the
#' clusters are unlikely to provide meaningful results on which the user can
#' better understand the inherent structure within the data.
#' 
#' @param x a matrix containing the data to cluster. The rows are the sample
#' observations, and the columns are the features.
#' @param K the number of clusters
#' @param prob a vector of probabilities to generate each cluster label. If
#' NULL, each cluster label has an equal chance of being selected.
#' @return a vector of clustering labels for each observation in 'x'.
random_clustering <- function(x, K, prob = NULL) {
  if (!is.null(prob)) {
    if (!is.numeric(prob)) {
      stop("The vector 'prob' must be 'numeric'.")
    }
    if (K != length(prob)) {
      stop("The length of 'prob' must equal 'K'.")
    }
    if (sum(prob) != 1) {
      stop("The sum of the probabilities must sum to 1.")
    }
    if (any(prob <= 0) || any(prob >= 1)) {
      stop("The cluster probabilties must be between 0 and 1.")
    }
  }

  sample(x = seq_len(K), size = nrow(x), replace = TRUE, prob = prob)
}
