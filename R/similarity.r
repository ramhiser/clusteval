#' A wrapper for the similarity measures
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @param method TODO
#' @return TODO
cluster_similarity <- function(cl1, cl2, method) {
	similarity_measures <- c("rand", "jaccard", "adjustedrand")
	method <- match.arg(method, similarity_measures)
	method <- paste(method, "_wrapper", sep = "")

  # TODO: For now, I have short-circuited this function
  # and only use the Jaccard score because it is the only one
  # of interest in our current paper.
  jaccard(cl1, cl2)	
}

#' Test if all the elements of the vector x are equal.
#'
#' TODO
#'
#' My code is based on the following Stack Overflow post.
#' http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
#' 
#' @param x TODO
#' @param tol TODO
#' @return TODO
vec_equal <- function(x,  tol = .Machine$double.eps ^ 0.5) {
  diff(range(x)) < tol
}

#' Co-Membership
#'
#' For a set of clustering labels, computes the co-membership of all
#' pairs of observations. Basically, the co-membership is defined
#' as the pairs of observations that are clustered together.
#'
#' TODO: Add purpose. (Label switching problem).
#' TODO: Add Kapp and Tibshirani (2007) Biostatistics paper as reference on 'co-membership' term.
#' TODO: Add to the name of the function.
#'
#' @param clust_labels a vector of clustering labels
#' @return a vector of co-membership bits
#' @examples
#' TODO
comembership <- function(clust_labels) {
  as.integer(combn(clust_labels, 2, vec_equal))
}

#' Summary of pairs of observations from two different clusterings of the same data set.
#'
#' TODO
#'
#' We assume that the lengths of each vector of cluster labels are equal.
#'
#' @export
#' @param cl1 the cluster labels
#' @param cl2 TODO
#' @return TODO
cluster_pairs <- function(cl1, cl2) {
  if (length(cl1) != length(cl2)) {
    stop("The length of 'cl1' must be equal to the length of 'cl2'.")
  }
  out <- list()
  out$pairs1 <- comembership(cl1)
  out$pairs2 <- comembership(cl2)
  out$and <- with(out, pairs1 & pairs2)
  out$or <- with(out, pairs1 | pairs2)
  out
}

#' Computes the Jaccard index for comparing two clusterings of the same data set.
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @Param cl2 TODO
#' @return numeric the Jaccard index for the two sets of cluster labels. If an error is encountered, we return NULL.
jaccard <- function(cl1, cl2) {
  cl_pairs <- cluster_pairs(cl1, cl2)
  jaccard_out <- plyr:::try_default(sum(cl_pairs$and) / (sum(cl_pairs$or)), default = NULL)

  # In the case where 'cl1' and 'cl2' contain all singletons, the Jaccard coefficient
  # results in the expression 0 / 0, which yields a NaN value in R.
  # We define such cases as 0.
  if (is.nan(jaccard_out)) {
    jaccard_out <- 0
  }
  jaccard_out
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

#' Wrapper for the Adjusted Rand index for comparing two clusterings of the same data set.
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @return TODO
adjustedrand_wrapper <- function(cl1, cl2) {
	adjustedRand(cl1, cl2, randMethod = "HA") # from 'clues' package on CRAN
}
