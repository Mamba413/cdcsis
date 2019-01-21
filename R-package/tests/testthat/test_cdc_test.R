library(testthat)
library(cdcsis)


test_that("Test Conditional-Test (Type-I error)", {
  num <- 50
  
  set.seed(123)
  cov_mat <- matrix(c(1, 0.36, 0.6, 0.36, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
  dat <- mvtnorm::rmvnorm(n = num, sigma = cov_mat)
  x <- dat[, 1]
  y <- dat[, 2]
  z <- dat[, 3]
  expect_gt(cdcov.test(x, y, z)[["p.value"]], 0.05)

  set.seed(123)
  z <- rnorm(num)
  x <- 0.5 * (z^3 / 7 + z / 2) + tanh(rnorm(num))
  x <- x + x^3 / 3
  y <- (z^3 + z) / 3 + rnorm(num)
  y <- y + tanh(y / 3)
  expect_gt(cdcov.test(x, y, z)[["p.value"]], 0.05)
})


test_that("Test Conditional-Test (Type-II error)", {
  num <- 50
  
  set.seed(1)
  cov_mat <- matrix(c(1, 0.7, 0.6, 0.7, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
  dat <- mvtnorm::rmvnorm(n = num, sigma = cov_mat)
  x <- dat[, 1]
  y <- dat[, 2]
  z <- dat[, 3]
  expect_lt(cdcov.test(x, y, z, width = 0.5)[["p.value"]], 0.05)
  
  set.seed(1)
  z <- matrix(rt(num * 4, df = 2), nrow = num)
  x <- z
  y <- cbind(sin(z[, 1]) + cos(z[, 2]) + (z[, 3])^2 + (z[, 4])^2, 
             (z[, 1])^2 + (z[, 2])^2 + z[, 3] + z[, 4])
  z <- z[, 1:2]
  expect_lt(cdcov.test(x, y, z, seed = 2)[["p.value"]], 0.05)
})
