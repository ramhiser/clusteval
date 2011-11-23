\name{landon_boot}
\alias{landon_boot}
\title{The worker function for Landon's Alternative Method}
\usage{
  landon_boot(x, idx, K, cluster_method, similarity_method,
    with_replacement = TRUE, ...)
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
  distribution of Landon's alternative method.
}
\details{
  TODO
}

