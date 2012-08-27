\name{sim_data}
\alias{sim_data}
\title{Wrapper function to generate data from a variety of distributions.}
\usage{
  sim_data(family = c("uniform", "normal", "student"), ...)
}
\arguments{
  \item{family}{the family of distributions from which to
  generate data}

  \item{...}{optional arguments that are passed to the
  data-generating function}
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

