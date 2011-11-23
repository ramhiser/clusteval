\name{clustomit_parboot}
\alias{clustomit_parboot}
\title{A Parametric Bootstrap Approach to Evaluate Clustering Stability via Cluster Omission.}
\usage{
  clustomit_parboot(x, num_clusters, cluster_method,
    similarity_method, noise = 0.1, B = 100, ...)
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
  The cluster method should take the data as 'x' and the
  'num_clusters'.
}

