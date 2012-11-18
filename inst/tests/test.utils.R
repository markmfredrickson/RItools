################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("Utility Functions")

test_that("xbal formula method", {
  # create a quick xBalance result
  df <- data.frame(Z = rep(c(1,0), 10),
                   X = rnorm(20))

  xb <- xBalance(Z ~ X, data = df)

  expect_identical(Z ~ X, as.formula(xb))
})

