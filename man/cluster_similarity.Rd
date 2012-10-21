\name{cluster_similarity}
\alias{cluster_similarity}
\title{Computes the similarity between two clusterings of the same data set.}
\usage{
  cluster_similarity(labels1, labels2,
    similarity = c("jaccard", "rand"), ...)
}
\arguments{
  \item{labels1}{a vector of \code{n} clustering labels}

  \item{labels2}{a vector of \code{n} clustering labels}

  \item{similarity}{the similarity statistic to calculate}

  \item{...}{Additional arguments passed to the similarity
  function.}
}
\value{
  the similarity between the two clusterings
}
\description{
  This function is a wrapper around the similarity indices
  implemented in the \code{clusteval} package. Currently,
  we have implemented the following similarity statistics:
  \itemize{ \item Rand index \item Jaccard coefficient }
}
\details{
  Note that the computed similarity indices may be
  substantially different than other packages, including
  the \code{fpc} and \code{clv} packages on CRAN. By
  default, we compute improved estimators, which are a work
  in progress. To calculate the commonly employed
  estimators, use the argument \code{estimation =
  "standard"}.
}
\examples{
# Notice that the number of comemberships is 'n choose 2'.
iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
cluster_similarity(iris_kmeans, iris_hclust)
cluster_similarity(iris_kmeans, iris_hclust, estimation = "standard")
cluster_similarity(iris_kmeans, iris_hclust, similarity = "rand")
}

