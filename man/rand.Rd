\name{rand}
\alias{rand}
\title{Computes the Rand similarity coefficient of two vectors of cluster labels.}
\usage{
  rand(labels1, labels2,
    estimation = c("glmm", "standard"))
}
\arguments{
  \item{labels1}{a vector of \code{n} clustering labels}

  \item{labels2}{a vector of \code{n} clustering labels}

  \item{estimation}{the method used to estimate the Rand
  similarity coefficient. By default, we use a GLMM
  method.}
}
\value{
  the Rand index for the two sets of cluster labels
}
\description{
  TODO
}
\details{
  To compute the contingency table, we use the
  \code{\link{comembership_table}} function.
}
\examples{
# We generate K = 3 labels for each of n = 10 observations and compute the
# Rand similarity index between the two clusterings.
set.seed(42)
K <- 3
n <- 10
labels1 <- sample.int(K, n, replace = TRUE)
labels2 <- sample.int(K, n, replace = TRUE)
rand(labels1, labels2, estimation = "glmm")
rand(labels1, labels2, estimation = "standard")

# Here, we cluster the \\code{\\link{iris}} data set with the K-means and
# hierarchical algorithms using the true number of clusters, K = 3.
# Then, we compute the Rand similarity index between the two clusterings.
iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
rand(iris_kmeans, iris_hclust, estimation = "glmm")
rand(iris_kmeans, iris_hclust, estimation = "standard")
}

