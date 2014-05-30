context("Adjusted Rand Index")

# Test related to ramhiser/clusteval#26
test_that("adjusted_rand returns 1 when when both vectors are singleton clusters", {
  n <- 10
  labels1 <- rep(1, 10)
  labels2 <- rep(2, 10)
  adjrand_out <- adjusted_rand(labels1, labels2)
  expect_equal(adjrand_out, 1)
})
