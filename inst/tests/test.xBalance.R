################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("xBalance main function")

test_that("Groups", {

  # create a quick xBalance result
  df <- data.frame(Z = rep(c(1,0), 10),
                   X = rnorm(20),
                   Y = rnorm(20),
                   W = rnorm(20),
                   S = as.factor(rep(c("A", "B"), each = 10)))

  xb <- xBalance(Z ~ X + Y + W, data = df, strata = data.frame(factor("none"), df$S))

  expect_is(xb$groups, "array")
  expect_equal(dim(xb$groups), c(2, 3, 1))
  
  xb.grp <- xBalance(Z ~ X + Y + W, data = df, strata = data.frame(factor("none"), df$S), 
                     groups = list("XY" = c("X", "Y"), "YW" = c("Y", "W")))

  # groups is an array of strata by test by groups (including "all")
  expect_equal(dim(xb.grp$groups), c(2, 3, 3))

  # Error checking
  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("XW" = c("X", "W"))),
               "Unknown variable")

  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("XZ" = c("X", "Z"))),
               "Treatment")

  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("Empty" = c())),
               "Empty group")

  # more tests/features: variables that expand into dummies/etc, passing groups in fmla or letting items of `groups` be little forumlae

})

