\name{sim_student}
\alias{sim_student}
\title{Generates random variates from multivariate Student's t populations.}
\usage{
  sim_student(n = rep(25, 5), p = 50, df = rep(6, 5),
    delta = 0, Sigma = diag(p), seed = NULL)
}
\arguments{
  \item{n}{a vector (of length M) of the sample sizes for
  each population}

  \item{p}{the dimension of the multivariate Student's t
  distributions}

  \item{df}{a vector (of length M) of the degrees of
  freedom for each population}

  \item{delta}{the fixed distance between each population
  and the origin}

  \item{Sigma}{the common covariance matrix}

  \item{seed}{seed for random number generation (If NULL,
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
  We generate \eqn{n_m} observations \eqn{(m = 1, \ldots,
  M)} from each of \eqn{M} multivariate Student's t
  distributions such that the Euclidean distance between
  each of the means and the origin is equal and scaled by
  \eqn{\Delta \ge 0}.
}
\details{
  Let \eqn{\Pi_m} denote the \eqn{m}th population with a
  \eqn{p}-dimensional multivariate Student's t
  distribution, \eqn{T_p(\mu_m, \Sigma_m, c_m)}, where
  \eqn{\mu_m} is the population location vector,
  \eqn{\Sigma_m} is the positive-definite covariance
  matrix, and \eqn{c_m} is the degrees of freedom.

  Let \eqn{e_m} be the \eqn{m}th standard basis vector
  (i.e., the \eqn{m}th element is 1 and the remaining
  values are 0). Then, we define \deqn{\mu_m = \Delta
  \sum_{j=1}^{p/M} e_{(p/M)(m-1) + j}.} Note that \code{p}
  must be divisible by \code{M}. By default, the first 10
  dimensions of \eqn{\mu_1} are set to \code{delta} with
  all remaining dimensions set to 0, the second 10
  dimensions of \eqn{\mu_2} are set to \code{delta} with
  all remaining dimensions set to 0, and so on.

  We use a common covariance matrix \eqn{\Sigma_m = \Sigma}
  for all populations.

  For small values of \eqn{c_m}, the tails are heavier,
  and, therefore, the average number of outlying
  observations is increased.

  By default, we let \eqn{M = 5}, \eqn{\Delta = 0},
  \eqn{\Sigma_m = I_p}, and \eqn{c_m = 6}, \eqn{m = 1,
  \ldots, M}, where \eqn{I_p} denotes the \eqn{p \times p}
  identity matrix. Furthermore, we generate 25 observations
  from each population by default.

  For \eqn{\Delta = 0} and \eqn{c_m = c}, \eqn{m = 1,
  \ldots, M}, the \eqn{M} populations are equal.
}
\examples{
data_generated <- sim_student(n = 10 * seq_len(5), seed = 42)
dim(data_generated$x)
table(data_generated$y)

data_generated2 <- sim_student(p = 10, delta = 2, df = rep(2, 5))
table(data_generated2$y)
sample_means <- with(data_generated2,
                     tapply(seq_along(y), y, function(i) {
                            colMeans(x[i,])
                     }))
(sample_means <- do.call(rbind, sample_means))
}

