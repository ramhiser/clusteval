\name{sim_normal}
\alias{sim_normal}
\title{Generates multivariate normal populations with intraclass covariance matrices.}
\usage{
  sim_normal(n = rep(25, 5), p = 50,
    rho = 0.1 * seq.int(1, 9, by = 2), delta = 0,
    epsilon = 1, sigma2 = 1, seed = NULL)
}
\arguments{
  \item{n}{a vector (of length M) of the sample sizes for
  each population}

  \item{p}{the dimension of the multivariate normal
  populations}

  \item{rho}{a vector (of length M) of the intraclass
  constants for each population}

  \item{delta}{the fixed distance between each population}

  \item{epsilon}{the amount of noise to add to the
  population means.}

  \item{sigma2}{the coefficient of each intraclass
  covariance matrix}

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
  We generate \eqn{n_m} observations \eqn({m = 1, \ldots,
  M}) from each of \eqn{M} multivariate normal
  distributions such that the Euclidean distance between
  the means of each of the populations is a fixed constant,
  \eqn{\Delta > 0}.
}
\details{
  By default, we consider \eqn{M = 5} populations.

  Let \eqn{\Pi_m} denote the \eqn{m}th population with a
  \eqn{p}-dimensional multivariate normal distribution,
  \eqn{N_p(\mu_m, \Sigma_m)} with mean vector \eqn{\mu_m}
  and covariance matrix \eqn{\Sigma_m}. Also, let \eqn{z_1,
  \ldots, z_M \sim N_p(0_p, \epsilon I_p)} be independently
  distributed vectors of noise for a specified
  \eqn{\epsilon > 0}, where \eqn{0_p} is a \eqn{p \times 1}
  vector of zeros. Also, let \eqn{e_m} be the \eqn{m}th
  standard basis vector (i.e. the \eqn{m}th element is 1
  and the remaining values are 0). Then, we define
  \eqn{\mu_m = \Delta e_m + z_m}. We introduce noise to the
  means to allow for more realistic simulated data set.

  By default, we let \eqn{\Delta = 0} and \eqn{\epsilon =
  0.1}.

  Also, we considered intraclass covariance (correlation)
  matrices such that \eqn{\Sigma_m = \sigma^2 (1 - \rho_m)
  J_p + \rho_m I_p}, where \eqn{-(p-1)^{-1} < \rho_m < 1},
  \eqn{I_p} is the \eqn{p \times p} identity matrix, and
  \eqn{J_p} denotes the \eqn{p \times p} matrix of ones.

  By default, we let \eqn{\sigma^2 = 1}.

  We generate \eqn{n_m} observations from population
  \eqn{\Pi_m}. By default, we generate 25 observations from
  each population.
}
\examples{
TODO
}

