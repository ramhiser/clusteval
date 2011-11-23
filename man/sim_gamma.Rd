\name{sim_gamma}
\alias{sim_gamma}
\title{Generates multivariate Gamma data}
\usage{
  sim_gamma(n = 25, delta = 0, seed = NULL)
}
\arguments{
  \item{n}{sample size of each population}

  \item{delta}{the fixed distance between each population}

  \item{seed}{Seed for random number generation. (If NULL,
  does not set seed)}
}
\value{
  data.frame. The 'Population' column denotes the
  population from which the observation in each row was
  generated. The remaining columns in each row contain the
  generated observation.
}
\description{
  TODO
}
\examples{
TODO
}

