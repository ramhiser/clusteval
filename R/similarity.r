#' A wrapper for the similarity measures
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @param method TODO
#' @return TODO
cluster_similarity <- function(cl1, cl2, similarity = c("jaccard", "rand"),
                               method = "independence") {
	similarity <- match.arg(similarity)
  method <- match.arg(method)

  # Currently, we effectively ignore the `method` argument and only use the
  # independence statistics.
  switch(similarity,
         jaccard = jaccard_indep(cl1, cl2),
         rand = rand_indep(cl1, cl2))
}

#' Cluster Comembership
#'
#' For a set of clustering labels, computes the comembership of all
#' pairs of observations. Basically, the comembership is defined
#' as the pairs of observations that are clustered together.
#'
#' TODO: Add purpose. (Label switching problem).
#' TODO: Add Kapp and Tibshirani (2007) Biostatistics paper as reference on 'co-membership' term.
#' TODO: Add to the name of the function.
#'
#' We use the \code{Rcpp} package to improve the runtime speed of
#' the \code{comembership} function.
#'
#' @export
#' @param labels a vector of clustering labels
#' @return a vector of co-membership bits
#' @examples
#' TODO
comembership <- function(labels) {
	.Call("rcpp_comembership", labels, PACKAGE = "clusteval")
}

#' Cluster Comembership
#'
#' For a set of clustering labels, computes the comembership of all
#' pairs of observations. Basically, the comembership is defined
#' as the pairs of observations that are clustered together.
#'
#' TODO: Add purpose. (Label switching problem).
#' TODO: Add Kapp and Tibshirani (2007) Biostatistics paper as reference on 'co-membership' term.
#' TODO: Add to the name of the function.
#'
#' We use the \code{Rcpp} package to improve the runtime speed of
#' the \code{comembership} function.
#'
#' @export
#' @param labels a vector of clustering labels
#' @return a vector of co-membership bits
#' @examples
#' TODO
comembership_table <- function(labels1, labels2) {
  if (length(labels1) != length(labels2)) {
    stop("The two vectors of cluster labels must be of equal length.");
  }

	.Call("rcpp_comembership_table", labels1, labels2, PACKAGE = "clusteval")
}

#' Computes the Jaccard index for comparing two clusterings of the same data set.
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @return numeric the Jaccard index for the two sets of cluster labels. If an error is encountered, we return NULL.
jaccard_indep <- function(cl1, cl2) {
  com_table <- comembership_table(cl1, cl2)
  jaccard_out <- with(com_table, AA / (AA + AD + DA))

  # In the case where 'cl1' and 'cl2' contain all singletons, the Jaccard
  # coefficient results in the expression 0 / 0, which yields a NaN value in R.
  # We define such cases as 0.
  if (is.nan(jaccard_out)) {
    warning("The two clusterings contain all singletons. Returning 0...")
    jaccard_out <- 0
  }
  jaccard_out
}

#' Computes the Rand index for comparing two clusterings of the same data set.
#'
#' TODO
#'
#' @export
#' @param cl1 TODO
#' @param cl2 TODO
#' @return numeric the Rand index for the two sets of cluster labels. If an error is encountered, we return NULL.
rand_indep <- function(cl1, cl2) {
  com_table <- comembership_table(cl1, cl2)
  with(com_table, (AA + DD) / (AA + AD + DA + DD))
}
