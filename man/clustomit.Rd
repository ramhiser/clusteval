\name{clustomit}
\alias{clustomit}
\title{Cluster Stability Evaluation via Cluster Omission}
\usage{
  clustomit(x, K, cluster_method,
    similarity_method = "jaccard", weighted_mean = TRUE,
    B = 100, with_replacement = FALSE,
    use_multicore = FALSE, ncpus = 1, ...)
}
\arguments{
  \item{x}{data matrix with N observations (rows) and p
  features (columns).}

  \item{num_clusters}{TODO}

  \item{K}{TODO}

  \item{similarity_method}{TODO}

  \item{weighted_mean}{logical Should the aggergate
  similarity score for each bootstrap replication be
  weighted by the number of observations in each of the
  observed clusters?}

  \item{B}{TODO}

  \item{use_multicore}{TODO}

  \item{ncpus}{TODO}

  \item{...}{TODO}
}
\value{
  list of scores by omitted cluster
}
\description{
  TODO: Add more thorough documentation
}
\details{
  This is a bootstrapping approach to evaluating the
  stability of a clustering algorithm through the
  clustering admissibility criteria in Fisher and Van Ness
  (1971). In particular, the focus here is on the cluster
  omission admissibility criterion.

  We remove each cluster at a time. Then, we bootstrap from
  the given data after removing the cluster, and use a
  measure of similarity between the the cluster labels of
  the original data set and the cluster labels of the
  bootstrapped sample.

  The cluster method should take the data as 'x' and the
  'num_clusters'.
}

