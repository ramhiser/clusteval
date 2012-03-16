library('testthat')
library('clusteval')
library('mvtnorm')


context("Cluster Omission Stability (COS) method works with artificial data.")

test_that("The COS method works for artificial data set #1 without error", {
  # Simulation 1 - Low-dimensional, well-separated populations with simple
  # covariance structure
  # n = 40
  # p = 2
  # M = 4 (# of populations)
  # For m = 1 to M
  #   n_m = 10
  #   mu_m = c(10m, 10m)
  #   Sigma_m = I_2
  set.seed(42)
  M <- 4
  p <- 2
  sample_sizes <- rep(10, M)
  means <- replicate(p, 10 * seq.int(M))
  Sigma <- diag(p)
  x <- lapply(seq_len(M), function(m) {
    rmvnorm(n = sample_sizes[m], mean = means[m,], sigma = diag(p))
  })
  x <- do.call(rbind, x)

  out <- clustomit(x = x, K = 3, cluster_method = "kmeans", B = 20)
  out <- clustomit(x = x, K = 4, cluster_method = "kmeans", B = 20)
  out <- clustomit(x = x, K = 5, cluster_method = "kmeans", B = 20)
  out <- clustomit(x = x, K = 6, cluster_method = "kmeans", B = 20)
})

test_that("The COS method works for artificial data set #2 without error", {
  # Simulation 2 - Low-dimensional with two extreme outyling singletons
  # n = 47 total observations
  # p = 2
  # M = 5 (# of populations)
  # Sigma_m = I_2 for all m = 1 to M
  # n_1 = n_2 = n_3 = 15
  # n_4 = n_5 = 1
  # mu_1 = (0,0)'
  # mu_2 = (2,2)'
  # mu_3 = (-2,-2)'
  # mu_4 = (20, 20)'
  # mu_5 = (-20, -20)

  set.seed(42)
  M <- 5
  p <- 2
  sample_sizes <- c(15, 15, 15, 1, 1)
  means <- matrix(c(0, 0, 2, 2, -2, 2, 20, 20, -20, -20), ncol = 2, byrow = TRUE)
  Sigma <- diag(p)
  x <- lapply(seq_len(M), function(m) {
    rmvnorm(n = sample_sizes[m], mean = means[m,], sigma = diag(p))
  })
  x <- do.call(rbind, x)

  # K = 3 is a problem case. Why?
  #trace("clustomit_boot", quote(if (any(is.na(omit_similarities))) { browser() }), at = 6, print = F)  
  out <- clustomit(x = x, K = 3, cluster_method = "kmeans", B = 20)
  out <- clustomit(x = x, K = 4, cluster_method = "kmeans", B = 20)
  out <- clustomit(x = x, K = 5, cluster_method = "kmeans", B = 20)  
  out <- clustomit(x = x, K = 6, cluster_method = "kmeans", B = 20)
})
