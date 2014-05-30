#' Generates random variates from multivariate Student's t populations.
#'
#' We generate \eqn{n_k} observations \eqn{(k = 1, \ldots, K_0)} from each of
#' \eqn{K_0} multivariate Student's t distributions such that the Euclidean
#' distance between each of the means and the origin is equal and scaled by
#' \eqn{\Delta \ge 0}.
#'
#' Let \eqn{\Pi_k} denote the \eqn{k}th population with a \eqn{p}-dimensional
#' multivariate Student's t distribution, \eqn{T_p(\mu_k, \Sigma_k, c_k)}, where
#' \eqn{\mu_k} is the population location vector, \eqn{\Sigma_k} is the
#' positive-definite covariance matrix, and \eqn{c_k} is the degrees of freedom.
#'
#' Let \eqn{e_k} be the \eqn{m}th standard basis vector (i.e., the \eqn{k}th
#' element is 1 and the remaining values are 0). Then, we define \deqn{\mu_k =
#' \Delta \sum_{j=1}^{p/K_0} e_{(p/K_0)(k-1) + j}.} Note that \eqn{p} must be
#' divisible by \eqn{K_0}. By default, the first 10 dimensions of \eqn{\mu_1}
#' are set to \eqn{\Delta} with all remaining dimensions set to 0, the second 10
#' dimensions of \eqn{\mu_2} are set to \eqn{\Delta} with all remaining
#' dimensions set to 0, and so on.
#'
#' We use a common covariance matrix \eqn{\Sigma_k = \Sigma} for all populations.
#'
#' For small values of \eqn{c_k}, the tails are heavier, and, therefore, the
#' average number of outlying observations is increased.
#'
#' By default, we let \eqn{K_0 = 5}, \eqn{\Delta = 0}, \eqn{\Sigma_k = I_p}, and
#' \eqn{c_k = 6}, \eqn{k = 1, \ldots, K_0}, where \eqn{I_p} denotes the
#' \eqn{p \times p} identity matrix. Furthermore, we generate 25 observations
#' from each population by default.
#'
#' For \eqn{\Delta = 0} and \eqn{c_k = c}, \eqn{k = 1, \ldots, K_0}, the \eqn{K_0}
#' populations are equal.
#'
#' @param n a vector (of length M) of the sample sizes for each population
#' @param p the dimension of the multivariate Student's t distributions
#' @param df a vector (of length M) of the degrees of freedom for each population
#' @param delta the fixed distance between each population and the origin
#' @param Sigma the common covariance matrix
#' @param seed seed for random number generation (If NULL, does not set seed)
#' @return named list containing:
#' \describe{
#'   \item{x:}{A matrix whose rows are the observations generated and whose
#'   columns are the \code{p} features (variables)}
#'   \item{y:}{A vector denoting the population from which the observation in
#'   each row was generated.}
#' }
#' @export
#' @examples
#' data_generated <- sim_student(n = 10 * seq_len(5), seed = 42)
#' dim(data_generated$x)
#' table(data_generated$y)
#'
#' data_generated2 <- sim_student(p = 10, delta = 2, df = rep(2, 5))
#' table(data_generated2$y)
#' sample_means <- with(data_generated2,
#'                      tapply(seq_along(y), y, function(i) {
#'                             colMeans(x[i,])
#'                      }))
#' (sample_means <- do.call(rbind, sample_means))
sim_student <- function(n = rep(25, 5), p = 50, df = rep(6, 5), delta = 0,
                        Sigma = diag(p), seed = NULL) {
  # The number of populations
  K_0 <- length(n)

  if (delta < 0) {
    stop("The value for 'delta' must be a nonnegative constant.")
  }
  if (length(n) != length(df)) {
    stop("The length of the vectors 'n' and 'df' must be equal.")
  }
  if (p %% K_0 != 0) {
    stop("We require that 'p' be divisible by 'K_0'")
  }
  if(!is.null(seed)) {
    set.seed(seed)
  }

  # A matrix whose rows are the population centroids.
  centroids <- lapply(seq.int(K_0), function(k) {
    mu_k <- matrix(0, nrow = K_0, ncol = p / K_0)
    mu_k[k, ] <- 1
    mu_k
  })
  centroids <- delta * do.call(cbind, centroids)

  # Generates the data in a list of length K_0.
  # Then, we rbind the data together.
  x <- lapply(seq_len(K_0), function(k) {
    rmvt(n[k], sigma = Sigma, df = df[k], delta = centroids[k,])
  })
  x <- do.call(rbind, x)
  y <- do.call(c, sapply(seq_len(K_0), function(k) {
    rep.int(k, n[k])
  }, simplify = FALSE))
  list(x = x, y = y)
}

