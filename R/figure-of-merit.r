#' Figure of Merit Method for Estimating the Predictive Power of a Clustering
#' Algorithm
#'
#' We provide an implementation of the Figure of Merit (FOM) stability method
#' proposed by Yeung et al. (2001) for estimating the predictive power of a
#' clustering algorithm. The FOM statistic aggregates the average distance of
#' each observation to its cluster centroid (typically, the cluster mean) after
#' removing each observation in sequence, similar to the jackknife and
#' leave-one-out cross-validation methods. The aggregate FOM score is a function
#' of the average of the aggregated distances corresponding to each sample being
#' removed. For observed FOM scores for clusterings obtained from two distinct
#' clustering algorithms, the smaller score indicates that the corresponding
#' clustering algorithm is preferred and has more predictive power. However,
#' Yeung et al. (2001) note that the FOM statistic can be used only for relative
#' comparisons of clustering algorithms on the same data set for a specified
#' value of \code{K}.
#' 
#' We require a clustering algorithm function to be specified in the argument
#' \code{cluster_method}. The function given should accept at least two
#' arguments:
#' \describe{
#'   \item{x}{matrix of observations to cluster}
#'   \item{K}{the number of clusters to find}
#'   \item{...}{additional arguments that can be passed on}
#' }
#' Also, the function given should return only clustering labels for each
#' observation in the matrix \code{x}. The additional arguments specified in
#' \code{...} are useful if a wrapper function is used: see the example below for
#' an illustration.
#'
#' @export
#' @param x data matrix with \code{n} observations (rows) and \code{p} features
#' (columns)
#' @param K the number of clusters to find with the clustering algorithm
#' specified in \code{cluster_method}
#' @param cluster_method a character string or a function specifying the
#' clustering algorithm that will be used. The method specified is matched with
#' the \code{\link{match.fun}} function. The function given should return only
#' clustering labels for each observation in the matrix \code{x}.
#' @param adjusted If specified, the adjusted FOM is calculated.
#' @param ... additional arguments passed to the function specified in
#' \code{cluster_method}
#' @return object of class \code{figure_of_merit}, which contains a named list
#' with elements
#' \describe{
#'   \item{scores:}{vector of length \code{n} containing the individual FOM
#' scores for each observation removed in sequence}
#'   \item{aggregate:}{the aggregate FOM score. This value is \code{adjusted} if
#' specified}
#' }
#' @references Yeung K., Haynor D., and Ruzzo W. (2001), Validating Clustering
#' for Gene Expression Data, _Bioinformatics_, 17, 4, 309–318.
#' \url{http://bioinformatics.oxfordjournals.org/content/17/4/309.abstract}
#' @examples
#' # First, we create a wrapper function for the K-means clustering algorithm
#' # that returns only the clustering labels for each observation (row) in
#' # \code{x}.
#' kmeans_wrapper <- function(x, K, num_starts = 10, ...) {
#'   kmeans(x = x, centers = K, nstart = num_starts, ...)$cluster
#' }
#'
#' # For this example, we generate five multivariate normal populations with the
#' # \code{sim_data} function.
#' set.seed(42)
#' x <- sim_data("normal", delta = 1.5)$x
#'
#' fom_out <- figure_of_merit(x = x, K = 4, cluster_method = "kmeans_wrapper")
#' fom_out2 <- figure_of_merit(x = x, K = 5, cluster_method = kmeans_wrapper)
figure_of_merit <- function(x, K, cluster_method, adjusted = TRUE, ...) {
  x <- as.matrix(x)
  K <- as.integer(K)
  cluster_method <- match.fun(cluster_method)

  N <- nrow(x)
  fom_scores <- numeric(N)
  for (i in seq_len(N)) {
    labels <- cluster_method(x[-i, ], K = K, ...)
    distances <- tapply(seq_along(labels), labels, function(cluster_idx) {
      x_cluster <- x[-i, ][cluster_idx, ]
      if (is.vector(x_cluster)) {
        x_cluster <- matrix(x_cluster, nrow = 1)
      }
      dist2xbar <- as.matrix(dist(rbind(x_cluster, colMeans(x_cluster))))
      sum(tail(dist2xbar, n = 1)^2)
    })
    fom_scores[i] <- sqrt(sum(distances) / N)
  }
  aggregate_fom <- sum(fom_scores)
  if (adjusted) {
    aggregate_fom <- sqrt(N / (N - K)) * aggregate_fom
  }
  obj <- list(
    scores = fom_scores,
    aggregate = aggregate_fom
  )
  class(obj) <- "figure_of_merit"
  obj
}
