################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("Plotting Functions")

test_that("Basic plot", {
  # create a quick xBalance result
  set.seed(20121119)
  Z <- rep(c(0,1), 10)
  df <- data.frame(Z = Z,
                   X1 = rnorm(20, mean = Z*2),
                   X2 = rnorm(20, mean = Z*3),
                   X3 = rnorm(20, mean = Z * -1),
                   X4 = rnorm(20, mean = Z * -0.5))

  xb <- xBalance(Z ~ ., data = df)

  plot(xb)

  opts <- options(warn = 2)
  # has a an argument to make only a right-sided abs difference plot
  # it will be a warning if absolute isn't a parameter
  plot(xb, absolute = T)

  options(opts)

  # has an argument to order the variables (from bottom to bottom)
  plot(xb, which.vars = c("X1", "X2", "X4", "X3"))
 
  # but including other variables is an error
  expect_error(plot(xb, which.vars = c("X1", "X2", "X3", "X4", "X5")), "Unknown variable\\(s\\): X5")
})

