\name{consensus}
\alias{consensus}
\title{Consensus Clustering}
\usage{
  consensus(x, num_clusters,
    methods = c("hierarchical", "kmeans", "pam", "model", "diana"),
    resample = c("boot", "subsample", "perturb"),
    num_resamples = 100, subsample_size = 100,
    perturb_noise = 0.1, use_multicore = FALSE, ncpus = 1,
    ...)
}
\arguments{
  \item{x}{the data matrix}

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
  TODO TODO: Add to the list of clustering algorithms
}

