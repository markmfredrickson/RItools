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

test_that("filling in arrays", {
  f <- function(alpha, beta, gamma) { return(1000 * alpha + 100 * beta + gamma) }
  p <- list(alpha = 1:10, beta = 1:5, gamma = 1:3)
  
  res <- farray(f, p)

  expect_equal(dim(res), c(alpha = 10, beta = 5, gamma = 3))

  expect_equal(res[7,3,2], f(7, 3, 2))

  # should also work for a single arg

  res1 <- farray(function(i) { i^2 }, list(foo = c(1,2,3)))

  expect_is(res1, "array")

})
