\name{sim_normal}
\alias{sim_normal}
\title{Generates random variates from multivariate normal populations with intraclass
covariance matrices.}
\usage{
sim_normal(n = rep(25, 5), p = 50, rho = rep(0.9, 5), delta = 0,
  sigma2 = 1, seed = NULL)
}
\arguments{
  \item{n}{a vector (of length K_0) of the sample sizes for
  each population}

  \item{p}{the dimension of the multivariate normal
  populations}

  \item{rho}{a vector (of length K_0) of the intraclass
  constants for each population}

  \item{delta}{the fixed distance between each population
  and the origin}

  \item{sigma2}{the coefficient of each intraclass
  covariance matrix}

  \item{seed}{seed for random number generation (If
  \code{NULL}, does not set seed)}
}
\value{
named list containing: \describe{ \item{x:}{A matrix whose
rows are the observations generated and whose columns are
the \code{p} features (variables)} \item{y:}{A vector
denoting the population from which the observation in each
row was generated.} }
}
\description{
We generate \eqn{n_k} observations \eqn{(k = 1, \ldots,
K_0)} from each of \eqn{K_0} multivariate normal
distributions such that the Euclidean distance between each
of the means and the origin is equal and scaled by
\eqn{\Delta \ge 0}.
}
\details{
Let \eqn{\Pi_k} denote the \eqn{k}th population with a
\eqn{p}-dimensional multivariate normal distribution,
\eqn{N_p(\mu_k, \Sigma_k)} with mean vector \eqn{\mu_k} and
covariance matrix \eqn{\Sigma_k}. Also, let \eqn{e_k} be
the \eqn{k}th standard basis vector (i.e., the \eqn{k}th
element is 1 and the remaining values are 0). Then, we
define \deqn{\mu_k = \Delta \sum_{j=1}^{p/K_0}
e_{(p/K_0)(k-1) + j}.} Note that \eqn{p} must be divisible
by \eqn{K_0}. By default, the first 10 dimensions of
\eqn{\mu_1} are set to \eqn{\Delta} with all remaining
dimensions set to 0, the second 10 dimensions of
\eqn{\mu_2} are set to \eqn{\Delta} with all remaining
dimensions set to 0, and so on.

Also, we consider intraclass covariance (correlation)
matrices such that \eqn{\Sigma_k = \sigma^2 (1 - \rho_k)
J_p + \rho_k I_p}, where \eqn{-(p-1)^{-1} < \rho_k < 1},
\eqn{I_p} is the \eqn{p \times p} identity matrix, and
\eqn{J_p} denotes the \eqn{p \times p} matrix of ones.

By default, we let \eqn{K_0 = 5}, \eqn{\Delta = 0}, and
\eqn{\sigma^2 = 1}. Furthermore, we generate 25
observations from each population by default.

For \eqn{\Delta = 0} and \eqn{\rho_k = \rho}, \eqn{k = 1,
\ldots, K_0}, the \eqn{K_0} populations are equal.
}
\examples{
data_generated <- sim_normal(n = 10 * seq_len(5), seed = 42)
dim(data_generated$x)
table(data_generated$y)

data_generated2 <- sim_normal(p = 10, delta = 2, rho = rep(0.1, 5))
table(data_generated2$y)
sample_means <- with(data_generated2,
                     tapply(seq_along(y), y, function(i) {
                            colMeans(x[i,])
                     }))
(sample_means <- do.call(rbind, sample_means))
}

