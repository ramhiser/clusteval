context("Meila's (2007) Variation of Information")

test_that("VI is equal to values obtained manually", {
  K <- 3
  labels1 <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3)
  labels2 <- c(1, 2, 2, 2, 3, 3, 3, 3, 3, 3)
  n <- length(labels1)

  probs1 <- as.numeric(table(labels1) / n)
  probs2 <- as.numeric(table(labels2) / n)

  # Entropy
  H1 <- -drop(probs1 %*% log(probs1))
  H2 <- -drop(probs2 %*% log(probs2))
  
  # Mutual Information
  joint_probs <- matrix(table(labels1, labels2), nrow=K) / n
  mutual_information <- 0
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      log_term <- log(joint_probs[k1, k2] / probs1[k1] / probs2[k2])
      if (log_term == -Inf) {
        log_term <- 0
      }
      mutual_information <- mutual_information + joint_probs[k1, k2] * log_term
    }
  }

  expected_vi <- H1 + H2 - 2 * mutual_information

  expect_equal(variation_information(labels1, labels2),
               expected_vi)
})

test_that("VI satisfies the three metric axioms", {
  # 1. VI = 0 when the clusterings are identical
  set.seed(42)
  K <- 3
  n <- 30
  labels1 <- labels2 <- sample.int(K, n, replace=TRUE)
  expect_equal(variation_information(labels1, labels2),
               0)

  # 2. Symmetry
  # vi(labels1, labels2) == vi(labels2, labels2)
  set.seed(42)
  K <- 3
  n <- 30
  labels1 <- sample.int(K, n, replace=TRUE)
  labels2 <- sample.int(K, n, replace=TRUE)
  expect_equal(variation_information(labels1, labels2),
               variation_information(labels2, labels1))

  # 3. Triangle Inequality
  # vi(labels1, labels2) + vi(labels2, labels3) >=  vi(labels1, labels3)
  set.seed(42)
  K <- 3
  n <- 30
  labels1 <- sample.int(K, n, replace=TRUE)
  labels2 <- sample.int(K, n, replace=TRUE)
  labels3 <- sample.int(K, n, replace=TRUE)
  vi1 <- variation_information(labels1, labels2)
  vi2 <- variation_information(labels2, labels3)
  vi3 <- variation_information(labels1, labels3)
  
  expect_true(vi1 + vi2 >= vi3)
})

test_that("The upper bounds for VI hold", {
  # Equation (24): vi(labels1, labels2) <= log(n)
  set.seed(42)
  K <- 3
  n <- 100
  labels1 <- sample.int(K, n, replace=TRUE)
  labels2 <- sample.int(K, n, replace=TRUE)
  expect_true(variation_information(labels1, labels2) <= log(n))
})
