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
	similarity_measures <- c("rand", "jaccard", "adjustedrand")
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