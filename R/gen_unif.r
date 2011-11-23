#' Generates four bivariate uniform populations.
#'
#' We generate 'n' observations from each of four bivariate
#' distributions such that the Euclidean distance between
#' each of the populations is a fixed constant, 'delta' > 0.
#'
#' To define the populations, let \eqn{X_1 \sim U(a_1, b_1)}
#' and \eqn{X_2 \sim U(a_2, b_2)} be independently distributed
#' uniform random variables with \eqn{a_1 < b_1} and \eqn{a_2 < b_2}.
#' Then, we denote the (independent) bivariate uniform distribution
#' as \eqn{X = (X_1, X_2)' \sim U(a_1, b_1) \times U(a_2, b_2)}.
#' Let \eqn{Pi_m} denote the \eqn{m}th population \eqn{(m = 1, \ldots, 4)}.
#' Then, we have the four populations:
#' \deqn{
#'  \Pi_1 = U(-1/2, 1/2) \times U(\Delta - 1/2, \Delta + 1/2),
#'  \Pi_2 = U(\Delta - 1/2, \Delta + 1/2) \times U(-1/2, 1/2),
#'  \Pi_3 = U(-1/2, 1/2) \times U(-\Delta - 1/2, -\Delta + 1/2),
#'  \Pi_4 = U(-\Delta - 1/2, -\Delta + 1/2) \times U(-1/2, 1/2).
#' }
#'
#' @param n sample size of each population
#' @param delta the fixed distance between each population
#' @param seed Seed for random number generation. (If NULL, does not set seed)
#' @return data.frame. The 'Population' column denotes the population from which
#' the observation in each row was generated. The remaining columns in each row
#' contain the generated observation.
#' @export
#' @examples
#' x <- gen_unif(50, delta = 1.5)
#' plyr:::ddply(x, .(Population), summarize, xbar_1 = mean(x1), xbar_2 = mean(x2))
gen_unif <- function(n = 25, delta = 0, seed = NULL) {
  if (delta < 0) {
    stop("The value for 'delta' must be a nonnegative constant.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  pop1 <- c(-1/2, 1/2, delta - 1/2, delta + 1/2)
  pop2 <- c(delta - 1/2, delta + 1/2, -1/2, 1/2)
  pop3 <- c(-1/2, 1/2, -delta - 1/2, -delta + 1/2)
  pop4 <- c(-delta - 1/2, -delta + 1/2, -1/2, 1/2)
  
  unif_pops <- rbind.data.frame(pop1, pop2, pop3, pop4)
  colnames(unif_pops) <- c("a1", "b1", "a2", "b2")
  unif_pops$n <- n
  
  bivar_unif <- function(n, a1, b1, a2, b2) {
    cbind(runif(n, a1, b1), runif(n, a2, b2))
  }
  
  x <- mdply(unif_pops, as.data.frame(bivar_unif), .expand = F)
  colnames(x) <- c("Population", "x1", "x2")
  x$Population <- as.factor(x$Population)  
  x
}