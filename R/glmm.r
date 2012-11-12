#' Generates two sequences of comemberships from a GLMM fit.
#'
#' TODO
#'
#' @export
#' @param n the number of observations to generate
#' @param alpha the random effects intercept
#' @param beta the random effects "slope"
#' @param sigma the random effect standard deviation
#' @return list of two binary sequences -- each of length \code{n}
generate_comemberships <- function(n, alpha, beta, sigma) {
  # Generate random effect.
  gamma <- rnorm(n, 0, sigma)

  # Generates two artificial sets of comemberships.
  # Really, we generate correlated binomial sequences.
  list(
    seq1 = rbinom(n = n, size = 1,
      prob = exp(alpha + gamma) / (1 + exp(alpha + gamma))),
    seq2 = rbinom(n = n, size = 1,
      prob = exp(alpha + beta + gamma) / (1 + exp(alpha + beta + gamma)))
  )
}
