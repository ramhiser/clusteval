\name{clustomit_boot}
\alias{clustomit_boot}
\title{The worker function for Cluster Omission}
\usage{
  clustomit_boot(x, idx, K, cluster_method,
    similarity_method, with_replacement = FALSE,
    obs_clusters, ...)
}
\arguments{
  \item{x}{data matrix with N observations (rows) and p
  features (columns).}

  \item{idx}{TODO}

  \item{K}{TODO}

  \item{cluster_method}{TODO}

  \item{similarity_method}{TODO}

  \item{with_replacement}{TODO}

  \item{obs_clusters}{vector cluster labels for each
  observation from the observed (original) data set}

  \item{...}{TODO}
}
\value{
  vector the similarity scores for each omitted cluster
}
\description{
  This function is repeatedly called by boot() from the
  'boot' package to generate the approximate sampling
  distribution of the similarity index of the cluster
  omission method.
}
\details{
  TODO
}

