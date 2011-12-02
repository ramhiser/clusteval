\name{consensus_boot}
\alias{consensus_boot}
\title{The worker function for consensus clustering}
\usage{
  consensus_boot(x, idx, num_clusters,
    methods = c("hierarchical", "kmeans", "pam", "model", "diana"),
    ...)
}
\arguments{
  \item{x}{the data matrix}

  \item{idx}{the row indices used with the bootstrap
  package, 'boot'.}

  \item{num_clusters}{the number of clusters}

  \item{methods}{a vector of clustering methods to apply in
  consensus}

  \item{...}{(optional) additional arguments to pass to the
  clustering methods}
}
\value{
  TODO
}
\description{
  TODO
}

