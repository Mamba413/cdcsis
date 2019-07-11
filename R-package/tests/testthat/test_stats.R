library(testthat)
library(cdcsis)

test_that("Test CDC statistics (Correlation) value", {
  set.seed(1)
  x <- rnorm(25)
  y <- rnorm(25)
  z <- rnorm(25)
  target_value <- 0.1810813
  names(target_value) <- "cdcor"
  expect_equal(cdcor(x, y, z, width = 1.0)[["statistic"]], target_value)
})
