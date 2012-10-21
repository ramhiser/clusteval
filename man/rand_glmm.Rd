\name{rand_glmm}
\alias{rand_glmm}
\title{Computes the Rand similarity coefficient of two clusterings using a
generalized linear mixed models (GLMM) approach.}
\usage{
  rand_glmm(labels1, labels2, glmm_obj = NULL,
    num_reps = 10000)
}
\arguments{
  \item{labels1}{a vector of \code{n} clustering labels}

  \item{labels2}{a vector of \code{n} clustering labels}

  \item{glmm_obj}{an optional object of class \code{glmmML}
  from which the Rand similarity coefficient is computed.
  If specified, the cluster label arguments are ignored.}

  \item{num_reps}{the number of random effects sampled}
}
\value{
  the estimated Rand coefficient for the two sets of
  cluster labels
}
\description{
  The GLMM estimator improves the estimation of the Rand
  similarity coefficient and involves a random effects
  term.
}
\examples{
\dontrun{
# We generate K = 3 labels for each of n = 10 observations and compute the
# Rand similarity coefficient between the two clusterings.
set.seed(42)
K <- 3
n <- 10
labels1 <- sample.int(K, n, replace = TRUE)
labels2 <- sample.int(K, n, replace = TRUE)
rand_glmm(labels1, labels2)

# Here, we cluster the \\code{\\link{iris}} data set with the K-means and
# hierarchical algorithms using the true number of clusters, K = 3.
# Then, we compute the Rand similarity coefficient between the two
# clusterings.
iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
rand_glmm(iris_kmeans, iris_hclust)
}
}

