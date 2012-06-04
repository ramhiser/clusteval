#' Generates multivariate Student's t populations.
#'
#' TODO
#'
#' TODO: Unit test:
#'  all(sim_data("student", seed = 42) == sim_data("student", seed = 42))
#'
#' @param n a vector (of length M) of the sample sizes for each population
#' @param p the dimension of the distributions
#' @param df a vector (of length M) of the degrees of freedom for each population
#' @param delta the fixed distance between each population
#' @param Sigma the common covariance matrix
#' @param seed Seed for random number generation. (If NULL, does not set seed)
#' @return data.frame. The 'Population' column denotes the population from which
#' the observation in each row was generated. The remaining columns in each row
#' contain the generated observation.
#' @export
#' @examples
#' TODO
sim_student <- function(n = c(15, 15, 15, 30, 30), p = 10, df = c(10, 10, 10, 3, 3), delta = 0,
  Sigma = diag(p), seed = NULL) {
  if (delta < 0) {
    stop("The value for 'delta' must be a nonnegative constant.")
  }
  if (length(n) != length(df)) {
    stop("The length of the vectors 'n' and 'df' must be equal.")
  }
  if(!is.null(seed)) {
    set.seed(seed)
  }

  # The number of populations
  M <- length(n)

  # A matrix whose rows are the population centroids.
  centroids <- lapply(seq_len(M), function(m) {
    e_m <- rep(0, p)
    e_m[c(m, 2*m)] <- delta
    e_m
  })


  # Generates the data in a list of length M.
  # Then, we rbind the data together.
  x <- lapply(seq_len(M), function(m) cbind(m, rmvt(n[m], sigma = Sigma, df = df[m], delta = centroids[[m]])))
  x <- do.call(rbind.data.frame, x)
  
  colnames(x) <- c("Population", paste("x", seq.int(ncol(x) - 1), sep = ""))
  x$Population <- as.factor(x$Population)  
  x
}

