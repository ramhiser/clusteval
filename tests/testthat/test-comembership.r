context("Comembership")

test_that("comemberships are calculated correctly for integers", {
  set.seed(42)
  K <- 3
  n <- 5
  labels <- sample.int(K, n, replace = TRUE)
  comembership_out <- comembership(labels)

  expected_comemberships <- c(1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
  expect_equal(comembership_out, expected_comemberships)
})

test_that("comemberships are calculated correctly for factors", {
  set.seed(42)
  K <- 3
  n <- 5
  labels <- sample.int(K, n, replace = TRUE)
  labels <- factor(labels)
  comembership_out <- comembership(labels)

  expected_comemberships <- c(1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
  expect_equal(comembership_out, expected_comemberships)
})

# Test related to ramhiser/clusteval#35
test_that("comemberships are calculated correctly for characters", {
  set.seed(42)
  K <- 3
  n <- 5
  labels <- sample.int(K, n, replace = TRUE)
  labels <- as.character(labels)
  comembership_out <- comembership(labels)

  expected_comemberships <- c(1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
  expect_equal(comembership_out, expected_comemberships)
})
