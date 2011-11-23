\name{sim_gamma}
\alias{sim_gamma}
\title{Generates four bivariate uniform populations.}
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
  We generate 'n' observations from each of four bivariate
  distributions such that the Euclidean distance between
  each of the populations is a fixed constant, 'delta' > 0.
}
\details{
  To define the populations, let \eqn{X_1 \sim U(a_1, b_1)}
  and \eqn{X_2 \sim U(a_2, b_2)} be independently
  distributed uniform random variables with \eqn{a_1 < b_1}
  and \eqn{a_2 < b_2}. Then, we denote the (independent)
  bivariate uniform distribution as \eqn{X = (X_1, X_2)'
  \sim U(a_1, b_1) \times U(a_2, b_2)}. Let \eqn{Pi_m}
  denote the \eqn{m}th population \eqn{(m = 1, \ldots, 4)}.
  Then, we have the four populations: \deqn{ \Pi_1 =
  U(-1/2, 1/2) \times U(\Delta - 1/2, \Delta + 1/2), \Pi_2
  = U(\Delta - 1/2, \Delta + 1/2) \times U(-1/2, 1/2),
  \Pi_3 = U(-1/2, 1/2) \times U(-\Delta - 1/2, -\Delta +
  1/2), \Pi_4 = U(-\Delta - 1/2, -\Delta + 1/2) \times
  U(-1/2, 1/2). }
}
\examples{
TODO
}

