\name{sample_vector}
\alias{sample_vector}
\title{Wrapper function to sample from a vector without any convenience features}
\usage{
sample_vector(x, ...)
}
\arguments{
  \item{x}{vector}

  \item{...}{additional arguments passed to
  \code{\link[base]{sample}}}
}
\value{
a vector containing a random permutation of \code{x}
}
\description{
This function is a wrapper for the base
\code{\link[base]{sample}} function that circumvents the
_undesired behavior_ from the _convenience feature_
discussed in the function's documentation.
}
\details{
For example, if \code{x} contains a single element, say,
42, we wish \code{sample(42)} to return only \code{42}
rather than a random permutation of \code{1:42}.
}

