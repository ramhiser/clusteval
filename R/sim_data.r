#' Wrapper function to generate data from a variety of distributions.
#'
#' TODO
#'
#' @param family the family of distributions from which to generate data
#' @param ... optional arguments that are passed to the data-generating function
#' @return data.frame. The 'Population' column denotes the population from which
#' the observation in each row was generated. The remaining columns in each row
#' contain the generated observation.
#' @export
#' @examples
#' TODO
sim_data <- function(family = c("uniform", "normal", "student", "gamma"), ...) {
  switch(family,
    uniform = sim_unif(...),
    normal = sim_normal(...),
    student = sim_student(...),
    gamma = sim_gamma(...)
  )
}