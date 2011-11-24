\name{sim_student}
\alias{sim_student}
\title{Generates multivariate Student's t populations.}
\usage{
  sim_student(n = c(25, 25, 50, 50), p = 10,
    df = c(10, 10, 3, 3), delta = 0, Sigma = diag(p),
    seed = NULL)
}
\arguments{
  \item{n}{a vector (of length M) of the sample sizes for
  each population}

  \item{p}{the dimension of the distributions}

  \item{df}{a vector (of length M) of the degrees of
  freedom for each population}

  \item{delta}{the fixed distance between each population}

  \item{Sigma}{the common covariance matrix}

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
\details{
  TODO: Unit test: all(sim_data("student", seed = 42) ==
  sim_data("student", seed = 42))
}
\examples{
TODO
}

