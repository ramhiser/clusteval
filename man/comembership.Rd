\name{comembership}
\alias{comembership}
\title{Cluster Comembership}
\usage{
  comembership(labels)
}
\arguments{
  \item{labels}{a vector of clustering labels}
}
\value{
  a vector of co-membership bits
}
\description{
  For a set of clustering labels, computes the comembership
  of all pairs of observations. Basically, the comembership
  is defined as the pairs of observations that are
  clustered together.
}
\details{
  TODO: Add purpose. (Label switching problem). TODO: Add
  Kapp and Tibshirani (2007) Biostatistics paper as
  reference on 'co-membership' term. TODO: Add to the name
  of the function.

  We use the \code{Rcpp} package to improve the runtime
  speed of the \code{comembership} function.
}
\examples{
TODO
}

