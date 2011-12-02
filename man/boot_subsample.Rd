\name{boot_subsample}
\alias{boot_subsample}
\title{Randomly subsample a matrix or data frame.}
\usage{
  boot_subsample(x, subsample_size)
}
\arguments{
  \item{x}{the data matrix}

  \item{subsample_size}{the number of observations (rows)
  to select at random from x.}
}
\value{
  a random subsample of the data matrix, x.
}
\description{
  To randomly subsample, we are using the parametric option
  in the boot function within the boot package. The
  'parametric' option requires the specification of a
  'ran.gen' function to generate observations based on the
  original data and a list of maximum likelihood estimators
  (MLE). We utilize this method, but instead of the MLE, we
  instead pass the subsample_size.
}
\details{
  As noted in the boot documentation: Use of sim =
  "parametric" with a suitable ran.gen allows the user to
  implement any types of nonparametric resampling which are
  not supported directly.

  TODO: If 'x' and 'subsample_size' break, it may be
  necessary to use 'data' and 'mle' instead.
}
\examples{
TODO
}

