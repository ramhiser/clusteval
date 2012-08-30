\name{jaccard_indep}
\alias{jaccard_indep}
\title{Computes the Jaccard similarity coefficient of two clusterings of the same
data set under the assumption that the two clusterings are independent.}
\usage{
  jaccard_indep(labels1, labels2)
}
\arguments{
  \item{labels1}{a vector of \code{n} clustering labels}

  \item{labels2}{a vector of \code{n} clustering labels}
}
\value{
  the Jaccard coefficient for the two sets of cluster
  labels (See Details.)
}
\description{
  For two clusterings of the same data set, this function
  calculates the Jaccard similarity coefficient of the
  clusterings from the comemberships of the observations.
  Basically, the comembership is defined as the pairs of
  observations that are clustered together.
}
\details{
  To calculate the Rand index, we compute the 2x2
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

  The Jaccard similarity coefficient is defined as: \deqn{J
  = \frac{n_{11}}{n_{11} + n_{10} + n_{01}}}.

  In the special case that the Jaccard coefficient results
  in \eqn{0/0}, we define \eqn{J = 0}. For instance, this
  case can occur when both clusterings consist of all
  singleton clusters.

  To compute the contingency table, we use the
  \code{\link{comembership_table}} function.
}
\examples{
\dontrun{
# We generate K = 3 labels for each of n = 10 observations and compute the
# Jaccard similarity coefficient between the two clusterings.
set.seed(42)
K <- 3
n <- 10
labels1 <- sample.int(K, n, replace = TRUE)
labels2 <- sample.int(K, n, replace = TRUE)
jaccard_indep(labels1, labels2)

# Here, we cluster the \\code{\\link{iris}} data set with the K-means and
# hierarchical algorithms using the true number of clusters, K = 3.
# Then, we compute the Jaccard similarity coefficient between the two
# clusterings.
iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
jaccard_indep(iris_kmeans, iris_hclust)
}
}

