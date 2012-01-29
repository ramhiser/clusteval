\name{cluster}
\alias{cluster}
\title{Cluster Data with a Selected Clustering Algorithm}
\usage{
  cluster(x, K, method = "kmeans", ...)
}
\arguments{
  \item{x}{a matrix containing the data to cluster. The
  rows are the sample observations, and the columns are the
  features.}

  \item{K}{the number of clusters}

  \item{method}{a string specifying the clustering
  algorithm to apply. Possible values are: 'diana',
  'kmeans', 'mclust', 'pam', and 'random'. See the Details
  section below.}

  \item{...}{optional arguments passed to the clustering
  algorithm}
}
\value{
  an object of class 'cluster_results'
}
\description{
  This function is a wrapper for a wide variety of
  clustering algorithms. For a specified clustering
  algorithm, we apply the algorithm to a given data set for
  a specified number of clusters
}
\details{
  The 'random' clustering algorithm randomly chooses one of
  K labels to apply to each observation. By default, the
  probabilities for each cluster are equal, but these can
  be customized.

  Custom algorithms are not yet allowed, but this feature
  is planned for a future release.
}

