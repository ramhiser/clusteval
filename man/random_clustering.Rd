\name{random_clustering}
\alias{random_clustering}
\title{Randomly cluster a data set into K clusters.}
\usage{
  random_clustering(x, K, prob = NULL)
}
\arguments{
  \item{x}{a matrix containing the data to cluster. The
  rows are the sample observations, and the columns are the
  features.}

  \item{K}{the number of clusters}

  \item{prob}{a vector of probabilities to generate each
  cluster label. If NULL, each cluster label has an equal
  chance of being selected.}
}
\value{
  a vector of clustering labels for each observation in
  'x'.
}
\description{
  For each observation (row) in 'x', one of K labels is
  randomly generated. By default, the probabilities of
  selecting each clustering label are equal, but this can
  be altered by specifying 'prob', a vector of
  probabilities for each cluster.
}
\details{
  Random clustering is often utilized as a baseline
  comparison clustering against which other clustering
  algorithms are employed to identify structure within the
  data. Typically, comparisons are made in terms of
  proposed clustering assessment and evaluation methods as
  well as clustering similarity measures. For the former, a
  specified clustering evaluation method is computed for
  the considered clustering algorithms as well as random
  clustering. If the clusters determined by a considered
  clustering algorithm do not differ significantly from the
  random clustering, we might conclude that the found
  clusters are no better than naively choosing clustering
  labels for each observation at random. Likewise, a
  similarity measure can be computed to compare the
  clusterings from each of a considered clustering
  algorithm and a random clustering: if the clusterings are
  significantly similar, once again, we might conclude the
  clusters found via the considered clustering algorithm do
  not differ significantly from those found at random. In
  either case, the clusters are unlikely to provide
  meaningful results on which the user can better
  understand the inherent structure within the data.
}

