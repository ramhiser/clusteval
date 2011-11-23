#' Generates multivariate Gamma data
#'
#' TODO
#'
#' @param n sample size of each population
#' @param delta the fixed distance between each population
#' @param seed Seed for random number generation. (If NULL, does not set seed)
#' @return data.frame. The 'Population' column denotes the population from which
#' the observation in each row was generated. The remaining columns in each row
#' contain the generated observation.
#' @examples
#' TODO
sim_gamma <- function(n = 25, delta = 0, seed = NULL) {
  if (delta < 0) {
    stop("The value for 'delta' must be a nonnegative constant.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  warning("This function is stubbed and has not yet been implemented.")
  NULL
}