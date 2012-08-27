\name{comembership_table}
\alias{comembership_table}
\title{Cluster Comembership}
\usage{
  comembership_table(labels1, labels2)
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

