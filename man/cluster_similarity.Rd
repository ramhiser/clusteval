\name{cluster_similarity}
\alias{cluster_similarity}
\title{Computes the similarity between two clusterings of the same data set.}
\usage{
  cluster_similarity(labels1, labels2,
    similarity = c("jaccard", "rand"),
    method = "independence")
}
\arguments{
  \item{labels1}{a vector of \code{n} clustering labels}

  \item{labels2}{a vector of \code{n} clustering labels}

  \item{similarity}{the similarity statistic to calculate}

  \item{method}{the model under which the statistic was
  derived}
}
\value{
  the similarity between the two clusterings
}
\description{
  For two clusterings of the same data set, this function
  calculates the similarity statistic specified of the
  clusterings from the comemberships of the observations.
  Basically, the comembership is defined as the pairs of
  observations that are clustered together.
}
\details{
  To calculate the similarity, we compute the 2x2
  contingency table, consisting of the following four
  cells: \describe{ \item{n_11}{the number of observation
  pairs where both observations are comembers in both
  clusterings} \item{n_10}{the number of observation pairs
  where the observations are comembers in the first
  clustering but not the second} \item{n_01}{the number of
  observation pairs where the observations are comembers in
  the second clustering but not the first} \item{n_00}{the
  number of observation pairs where neither pair are
  comembers in either clustering} }

  Currently, we have implemented the following similarity
  statistics: \itemize{ \item Rand index \item Jaccard
  coefficient }

  To compute the contingency table, we use the
  \code{\link{comembership_table}} function.
}
\examples{
# Notice that the number of comemberships is 'n choose 2'.
iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
cluster_similarity(iris_kmeans, iris_hclust)
}

