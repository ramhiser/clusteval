#' Computes the Rand similarity coefficient of two vectors of cluster labels.
#'
#' TODO
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @param estimation the method used to estimate the Rand similarity
#' coefficient. By default, we use a GLMM method.
#' @return the Rand index for the two sets of cluster labels
#' @export
#' @examples
#'
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Rand similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' rand(labels1, labels2, estimation = "glmm")
#' rand(labels1, labels2, estimation = "standard")
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Rand similarity index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' rand(iris_kmeans, iris_hclust, estimation = "glmm")
#' rand(iris_kmeans, iris_hclust, estimation = "standard")
#'
rand <- function(labels1, labels2, estimation = c("glmm", "standard")) {
  estimation <- match.arg(estimation)
    
  switch(estimation,
         glmm = rand_glmm(labels1, labels2),
         standard = rand_standard(labels1, labels2))
}

#' Computes the Rand similarity coefficient of two clusterings using a
#' generalized linear mixed models (GLMM) approach.
#'
#' The GLMM estimator improves the estimation of the Rand similarity
#' coefficient and involves a random effects term.
#'
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @param glmm_obj an optional object of class \code{glmmML} from which the
#' Rand similarity coefficient is computed. If specified, the cluster label
#' arguments are ignored.
#' @param num_reps the number of random effects sampled
#' @return the estimated Rand coefficient for the two sets of cluster labels
#' @export
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Rand similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' rand_glmm(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Rand similarity coefficient between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' rand_glmm(iris_kmeans, iris_hclust)
#' }
rand_glmm <- function(labels1, labels2, glmm_obj = NULL, num_reps = 10000) {
  warning("The GLMM estimator of the Rand index is experimental.")
  if (!is.null(glmm_obj)) {
    if (!inherits(glmm_obj, "glmmML")) {
      stop("The 'glmm_obj' must be of class 'glmmML'.")
    }
    alpha <- as.numeric(glmm_obj$coefficients[1])
    beta <- as.numeric(glmm_obj$coefficients[2])
    sigma <- glmm_obj$sigma
  } else {
    comem_glmm <- comembership_glmm(labels1 = labels1, labels2 = labels2)
    alpha <- comem_glmm$alpha
    beta <- comem_glmm$beta
    sigma <- comem_glmm$sigma
  }
  
  gamma <- rnorm(num_reps, 0, sd = sigma)
  p1 <- mean(exp(alpha + gamma) / (1 + exp(alpha + gamma)))
  p2 <- mean(exp(alpha + beta + gamma) / (1 + exp(alpha + beta + gamma)))
  p1 * p2 + (1 - p1) * (1 - p2)
}

#' Computes the standard Rand similarity coefficient of two clusterings of the
#' same data set.
#'
#' The standard estimator implicitly involves the following three assumptions:
#' \enumerate{
#'   \item TODO: Fixed prob
#'   \item TODO: Independence within
#'   \item TODO: Independence between
#' }
#'
#' For two clusterings of the same data set, this function calculates the Rand
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Rand index, we compute the 2x2 contingency table, consisting
#' of the following four cells:
#' \describe{
#'   \item{n_11}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' The Rand similarity index is defined as:
#' \deqn{R = \frac{n_{11} + n_{00}}{n_{11} + n_{10} + n_{01} + n_{00}}}.
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Rand coefficient for the two sets of cluster labels (See
#' Details.)
#' @export
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Rand similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' rand_standard(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Rand similarity coefficient between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' rand_standard(iris_kmeans, iris_hclust)
#' }
rand_standard <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, (n_11 + n_00) / (n_11 + n_10 + n_01 + n_00))
}
