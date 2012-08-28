\name{sim_unif}
\alias{sim_unif}
\title{Generates random variates from five multivariate uniform populations.}
\usage{
  sim_unif(n = rep(25, 5), delta = 0, seed = NULL)
}
\arguments{
  \item{n}{a vector (of length M = 5) of the sample sizes
  for each population}

  \item{delta}{the fixed distance between each population
  and the origin}

  \item{seed}{Seed for random number generation. (If NULL,
  does not set seed)}
}
\value{
  named list containing: \describe{ \item{x:}{A matrix
  whose rows are the observations generated and whose
  columns are the \code{p} features (variables)}
  \item{y:}{A vector denoting the population from which the
  observation in each row was generated.} }
}
\description{
  We generate \code{n} observations from each of four
  trivariate distributions such that the Euclidean distance
  between each of the populations is a fixed constant,
  \code{delta} > 0.
}
\details{
  To define the populations, let \eqn{x = (X_1, \ldots,
  X_p)'} be a multivariate uniformly distributed random
  vector such that \eqn{X_j \sim U(a_j, b_j)} is an
  independently distributed uniform random variable with
  \eqn{a_j < b_j} for \eqn{j = 1, \ldots, p}. Let
  \eqn{Pi_m} denote the \eqn{m}th population \eqn{(m = 1,
  \ldots, 5)}. Then, we have the five populations:
  \deqn{\Pi_1 = U(-1/2, 1/2) \times U(\Delta - 1/2, \Delta
  + 1/2) \times U(-1/2, 1/2) \times U(-1/2, 1/2),}
  \deqn{\Pi_2 = U(\Delta - 1/2, \Delta + 1/2) \times
  U(-1/2, 1/2) \times U(-1/2, 1/2) \times U(-1/2, 1/2),}
  \deqn{\Pi_3 = U(-1/2, 1/2) \times U(-\Delta - 1/2,
  -\Delta + 1/2) \times U(-1/2, 1/2) \times U(-1/2, 1/2),}
  \deqn{\Pi_4 = U(-1/2, 1/2) \times U(-1/2, 1/2) \times
  U(-\Delta - 1/2, -\Delta + 1/2) \times U(-1/2, 1/2),}
  \deqn{\Pi_5 = U(-1/2, 1/2) \times U(-1/2, 1/2) \times
  U(-1/2, 1/2) \times U(\Delta - 1/2, \Delta + 1/2).}

  We generate \eqn{n_m} observations from population
  \eqn{\Pi_m}.

  For \eqn{\Delta = 0} and \eqn{\rho_m = \rho}, \eqn{m = 1,
  \ldots, M}, the \eqn{M} populations are equal.

  Notice that the support of each population is a unit
  hypercube with 4 features. Moreover, for \eqn{\Delta \ge
  1}, the populations are mutually exclusive and entirely
  separated.
}
\examples{
data_generated <- sim_unif(seed = 42)
dim(data_generated$x)
table(data_generated$y)

data_generated2 <- sim_unif(n = 10 * seq_len(5), delta = 1.5)
table(data_generated2$y)
sample_means <- with(data_generated2,
                     tapply(seq_along(y), y, function(i) {
                            colMeans(x[i,])
                     }))
(sample_means <- do.call(rbind, sample_means))
}

