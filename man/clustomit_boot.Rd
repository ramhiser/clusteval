\name{clustomit_boot}
\alias{clustomit_boot}
\title{A Bootstrap Approach to Evaluate Clustering Stability via Cluster Omission.}
\usage{
  clustomit_boot(x, num_clusters, cluster_method,
    similarity_method, B = 100, ...)
}
\arguments{
  \item{x}{data matrix with N observations (rows) and p
  features (columns).}

  \item{num_clusters}{TODO}

  \item{cluster_method}{TODO}

  \item{similarity_method}{TODO}

  \item{noise}{TODO}

  \item{B}{TODO}

  \item{...}{TODO}
}
\value{
  list of scores by omitted cluster
}
\description{
  This is a bootstrapping approach to evaluating the
  stability of a clustering algorithm through the
  clustering admissibility criteria in Fisher and Van Ness
  (1971). In particular, the focus here is on the cluster
  omission admissibility criterion.
}
\details{
  We remove each cluster at a time. Then, we bootstrap from
  the given data after removing the cluster, and use a
  measure of similarity between the the cluster labels of
  the original data set and the cluster labels of the
  bootstrapped sample.

  The cluster method should take the data as 'x' and the
  'num_clusters'.
}

