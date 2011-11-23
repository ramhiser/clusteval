\name{ramey_boot}
\alias{ramey_boot}
\title{The worker function for Ramey's Alternative Method}
\usage{
  ramey_boot(x, idx, K, cluster_method, similarity_method,
    with_replacement = TRUE, clusters, ...)
}
\arguments{
  \item{x}{data matrix with N observations (rows) and p
  features (columns).}

  \item{idx}{TODO}

  \item{K}{TODO}

  \item{cluster_method}{TODO}

  \item{similarity_method}{TODO}

  \item{with_replacement}{TODO}

  \item{...}{TODO}
}
\value{
  list with results (TODO: Add more detail)
}
\description{
  This function is repeatedly called by boot() from the
  'boot' package to generate the approximate sampling
  distribution of Ramey's alternative method.
}
\details{
  TODO
}

