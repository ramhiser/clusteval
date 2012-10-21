\name{comembership_glmm}
\alias{comembership_glmm}
\title{For two vectors of clustering labels obtained from the same data set,
we fit a Generalized Linear Mixed Model (GLMM) to the comemberships
of each set of labels.}
\usage{
  comembership_glmm(labels1, labels2)
}
\arguments{
  \item{labels1}{a vector of \code{n} clustering labels}

  \item{labels2}{a vector of \code{n} clustering labels}
}
\value{
  list containing the fitted \code{glmmML} object, the
  comemberships, and the explicit GLMM parameter estimates
}
\description{
  We fit a GLMM to the comemberships using Gauss-Hermite
  quadrature via the \code{glmmML} function from the
  \code{glmmML} package.
}

