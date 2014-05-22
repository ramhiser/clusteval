#' Computes the similarity between two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the 
#' similarity statistic specified of the clusterings from the comemberships of
#' the observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the similarity, we compute the 2x2 contingency table, consisting
#' of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' Currently, we have implemented the following similarity statistics:
#' \itemize{
#'   \item Adjusted Rand index
#'   \item Dice coefficient
#'   \item Fowlkes-Mallows coefficient
#'   \item Jaccard coefficient
#'   \item Phi coefficient
#'   \item Rand index
#'   \item Rogers-Tanimoto coefficient
#'   \item Russel-Rao coefficient
#'   \item Sokal-Sneath coefficient
#' }
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @param similarity the similarity statistic to calculate. See
#' \code{\link{similarity_methods}} for a listing of available similarity
#' methods. By default, the adjusted Rand index is used.
#' @return the similarity between the two clusterings
#' @examples
#' # Notice that the number of comemberships is 'n choose 2'.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' cluster_similarity(iris_kmeans, iris_hclust)
cluster_similarity <- function(labels1, labels2, similarity = "adjusted_rand") {
	similarity <- match.arg(similarity, similarity_methods()$method)
  similarity_fun <- match.fun(similarity)
  similarity_fun(labels1, labels2)
}

#' Brief description of all similarity functions in the 'clusteval' package.
#'
#' For all similarity functions in the \code{\link{clusteval}} package, we
#' summarize them in a data frame.
#'
#' @export
#' @return data frame describing each similarity function in the package. The
#' columns of the data frame are:
#' \describe{
#'   \item{method}{the name of the R function to calculate the similarity method}
#'   \item{name}{The name of the similarity method}
#' }
similarity_methods <- function() {
  similarity_summary <- rbind.data.frame(
    c("adjusted_rand", "Adjusted Rand"),
    c("dice", "Dice"),
    c("fowlkes_mallows", "Fowlkes-Mallows"),
    c("jaccard", "Jaccard"),
    c("phi", "Phi"),
    c("rand", "Rand"),
    c("rogers_tanimoto", "Rogers-Tanimoto"),
    c("russel_rao", "Russel-Rao"),
    c("sokal_sneath", "Sokal-Sneath")
  )
  colnames(similarity_summary) <- c('method', 'name')
  similarity_summary$method <- as.character(similarity_summary$method)
  similarity_summary$name <- as.character(similarity_summary$name)

  similarity_summary
}


#' Computes the adjusted Rand similarity index of two clusterings of the same
#' data set.
#'
#' For two clusterings of the same data set, this function calculates the
#' adjusted Rand similarity coefficient of the clusterings from the
#' comemberships of the observations.
#'
#' The adjusted Rand index is a variant of the Rand index that is corrected for
#' chance. We refer the interested reader to the Wikipedia entry for an overview
#' of the formula:
#' \url{http://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index}
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the adjusted Rand index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # adjusted Rand index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' adjusted_rand(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the adjusted Rand index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' adjusted_rand(iris_kmeans, iris_hclust)
#' }
adjusted_rand <- function(labels1, labels2) {
  labels1 <- factor(as.vector(labels1))
  labels2 <- factor(as.vector(labels2))
  n <- length(labels1)
  if (n != length(labels2)) {
    stop("The two vectors of cluster labels must be of equal length.")
  }

  # Summarizes the contingency table of agreement
  table_out <- table(labels1, labels2)
  
  # When each vector is a singleton cluster, the Adjusted Rand should be 1.
  if (all(dim(table_out) == c(1, 1))) {
    return(1)
  }
  margin1 <- as.vector(margin.table(table_out, 1))
  margin2 <- as.vector(margin.table(table_out, 2))
  
  # Calculates the quantities employed in the adjusted Rand index
  margin1_sum <- sum(choose(margin1, 2))
  margin2_sum <- sum(choose(margin2, 2))
  index <- sum(choose(as.vector(table_out), 2))
  expected_index <- margin1_sum * margin2_sum / choose(n, 2)
  max_index <- (margin1_sum + margin2_sum) / 2

  (index - expected_index) / (max_index - expected_index)
}

#' Computes the Dice similarity index of two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the Dice
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Dice index, we compute the 2x2 contingency table, consisting
#' of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Dice similarity index is defined as:
#' \deqn{\frac{2 * n_{11}}{2 n_{11} + n_{10} + n_{01}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Dice index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Dice similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' dice(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Dice similarity index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' dice(iris_kmeans, iris_hclust)
#' }
dice <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, 2 * n_11 / (2 * n_11 + n_10 + n_01))
}

#' Computes the Fowlkes-Mallows similarity index of two clusterings of the same
#' data set.
#'
#' For two clusterings of the same data set, this function calculates the
#' Fowlkes-Mallows similarity coefficient of the clusterings from the
#' comemberships of the observations. Basically, the comembership is defined as
#' the pairs of observations that are clustered together.
#'
#' To calculate the Fowlkes-Mallows index, we compute the 2x2 contingency table, consisting
#' of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Fowlkes-Mallows similarity index is defined as:
#' \deqn{\frac{n_{11}}{\sqrt{(n_{11} + n_{10})(n_{11} + n_{01})}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Fowlkes-Mallows index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Fowlkes-Mallows similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' fowlkes_mallows(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Fowlkes-Mallows similarity index between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' fowlkes_mallows(iris_kmeans, iris_hclust)
#' }
fowlkes_mallows <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, n_11 / sqrt((n_11 + n_10) * (n_11 + n_01)))
}

#' Computes the Jaccard similarity coefficient of two clusterings of the same
#' data set.
#'
#' For two clusterings of the same data set, this function calculates the Jaccard
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Jaccard coefficient, we compute the 2x2 contingency table,
#' consisting of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Jaccard similarity coefficient is defined as:
#' \deqn{J = \frac{n_{11}}{n_{11} + n_{10} + n_{01}}.}
#'
#' In the special case that the Jaccard coefficient results in \eqn{0/0},
#' we define \eqn{J = 0}. For instance, this case can occur when both clusterings
#' consist of all singleton clusters.
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Jaccard coefficient for the two sets of cluster labels (See
#' Details.)
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Jaccard similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' jaccard(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Jaccard similarity coefficient between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' jaccard(iris_kmeans, iris_hclust)
#' }
jaccard <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  jaccard_out <- with(com_table, n_11 / (n_11 + n_10 + n_01))

  # In the case where 'labels1' and 'labels2' contain all singletons, the Jaccard
  # coefficient results in the expression 0 / 0, which yields a NaN value in R.
  # We define such cases as 0.
  if (is.nan(jaccard_out)) {
    warning("The two clusterings contain all singletons -- returning 0.")
    jaccard_out <- 0
  }
  jaccard_out
}

#' Computes the Phi coefficient of two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the Phi
#' coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Phi coefficient, we compute the 2x2 contingency table,
#' consisting of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Phi coefficient is defined as:
#' \deqn{\frac{n_{11} * n_{00} - n_{10} * n_{01}}{\sqrt{(n_{11} + n_{10})(n_{11} + n_{01})(n_{00} + n_{10})(n_{00} + n_{01})}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Phi index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Phi coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' phi(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Phi coefficient between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' phi(iris_kmeans, iris_hclust)
#' }
phi <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)

  numerator <- with(com_table, (n_11 * n_00) - (n_10 * n_01))
  denom <- with(com_table, sqrt((n_11 + n_10) * (n_11 + n_01) * (n_00 + n_10)
                                * (n_00 + n_01)))

  numerator / denom
}

#' Computes the Rand similarity index of two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the Rand
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Rand index, we compute the 2x2 contingency table, consisting
#' of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Rand similarity index is defined as:
#' \deqn{\frac{n_{11} + n_{00}}{n_{11} + n_{10} + n_{01} + n_{00}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Rand index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Rand similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' rand(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Rand similarity index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' rand(iris_kmeans, iris_hclust)
#' }
rand <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, (n_11 + n_00) / (n_11 + n_10 + n_01 + n_00))
}

#' Computes the Rogers-Tanimoto similarity of two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the
#' Rogers-Tanimoto similarity coefficient of the clusterings from the
#' comemberships of the observations. Basically, the comembership is defined as
#' the pairs of observations that are clustered together.
#'
#' To calculate the Rogers-Tanimoto similarity, we compute the 2x2 contingency table,
#' consisting of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Rogers-Tanimoto similarity is defined as:
#' \deqn{\frac{n_{11} + n_{00}}{n_{11} + 2 (n_{10} + n_{01}) + n_{00}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Rogers-Tanimoto index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Rogers-Tanimoto similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' rogers_tanimoto(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Rogers-Tanimoto similarity index between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' rogers_tanimoto(iris_kmeans, iris_hclust)
#' }
rogers_tanimoto <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, (n_11 + n_00) / (n_11 + 2 * (n_10 + n_01) + n_00))
}

#' Computes the Russel-Rao similarity of two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the Russel-Rao
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Russel-Rao similarity, we compute the 2x2 contingency table,
#' consisting of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Russel-Rao similarity is defined as:
#' \deqn{\frac{n_{11}}{n_{11} + n_{10} + n_{01} + n_{00}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Russel-Rao index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Russel-Rao similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' russel_rao(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Russel_Rao similarity index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' russel_rao(iris_kmeans, iris_hclust)
#' }
russel_rao <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, n_11 / (n_11 + n_10 + n_01 + n_00))
}

#' Computes the Sokal-Sneath similarity of two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the Sokal-Sneath
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Sokal-Sneath similarity, we compute the 2x2 contingency table,
#' consisting of the following four cells:
#' \describe{
#'   \item{n_11:}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10:}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01:}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00:}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Sokal-Sneath similarity is defined as:
#' \deqn{\frac{2 (n_{11} + n_{00})}{2 n_{11} + n_{10} + n_{01} + 2 n_{00}}.}
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Sokal-Sneath index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Sokal-Sneath similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' sokal_sneath(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Sokal_Sneath similarity index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' sokal_sneath(iris_kmeans, iris_hclust)
#' }
sokal_sneath <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, 2 * (n_11 + n_00) / (2 * n_11 + n_10 + n_01 + 2 * n_00))
}
