\name{cluster_wrapper_mult}
\alias{cluster_wrapper_mult}
\title{A Wrapper Function for a List of Clustering Algorithms}
\usage{
  cluster_wrapper_mult(x, K,
    methods = c("kmeans", "pam", "hierarchical"),
    parallel = FALSE, ...)
}
\arguments{
  \item{x}{the data matrix}

  \item{K}{the number of clusters}

  \item{methods}{a vector of clustering methods to apply in
  consensus}

  \item{parallel}{If TRUE, applies each clustering
  algorithm in parallel. See Details.}

  \item{...}{(optional) additional arguments to pass to the
  clustering methods}
}
\description{
  TODO TODO: Add to the list of clustering algorithms TODO:
  How to do everything in parallel? TODO: For now, to keep
  things simple, I am using the already made clustering
  algorithm As such, renamed function to
  cluster_wrapper_mult. Eventually, use one and only one.
  This should be straightforward, but it would break too
  much for now.
}

