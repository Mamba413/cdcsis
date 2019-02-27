library(testthat)
library(cdcsis)


test_that("Test Conditional-Screening (Univariate Explanation Variable)", {
  set.seed(1)
  num <- 100
  p <- 300
  x <- matrix(rnorm(num * p), nrow = num)
  z <- rnorm(num)
  y <- 3 * x[, 1] + 1.5 * x[, 2] + 4 * z * x[, 5] + rnorm(num)
  res <- cdcsis(x, y, z)
  expect_true(all(c(1, 2, 5) %in% head(res[["ix"]], n = 10)))
})


test_that("Test Conditional-Screening (Multivariate Explanation Variable)", {
  set.seed(1)
  num <- 100
  p <- 200
  x <- matrix(rnorm(num * p), nrow = num)
  z <- rnorm(num)
  y <- 3 * x[, 1] + 1.5 * x[, 2] + 4 * z * x[, 5] + rnorm(num)
  x <- as.list(as.data.frame(x))
  x <- lapply(x, as.matrix)
  x[[1]] <- cbind(x[[1]], x[[2]])
  x[[2]] <- NULL
  res <- cdcsis(x, y, z)
  expect_true(all(c(1, 4) %in% head(res[["ix"]], n = 10)))
})