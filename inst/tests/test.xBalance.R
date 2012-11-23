################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("xBalance main function")

test_that("Groups", {

  # create a quick xBalance result
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = rnorm(n),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   S = as.factor(rep(c("A", "B"), each = n/2)))

  xb <- xBalance(Z ~ X + Y + W, data = df, strata = data.frame(factor("none"), df$S))

  expect_is(xb$groups, "array")
  expect_equal(dim(xb$groups), c(2, 3, 1))
  expect_true(!any(is.na(xb$groups)))
  
  xb.grp <- xBalance(Z ~ X + Y + W, data = df, strata = data.frame(factor("none"), df$S), 
                     groups = list("XY" = c("X", "Y"), "YW" = c("Y", "W")))

  # groups is an array of strata by test by groups (including "all")
  expect_equal(dim(xb.grp$groups), c(2, 3, 3))
  expect_true(!any(is.na(xb.grp$groups)))

  # first off, the groups should have different values for the chisquared tests
  expect_true(!identical(xb.grp$groups[,"chisquare","All"], xb.grp$groups[,"chisquare","XY"]))
  # second, the degrees of freedom should be 3, 2, 2 respectively
  expect_equal(xb.grp$groups[1, "df", "All"], 3)
  expect_equal(xb.grp$groups[1, "df", "XY"], 2)
  expect_equal(xb.grp$groups[1, "df", "YW"], 2)

  # Error checking
  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("XW" = c("X", "W"))),
               "Unknown variable")

  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("XZ" = c("X", "Z"))),
               "Treatment")

  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("Empty" = c())),
               "Empty group")

  # more tests/features: variables that expand into dummies/etc, passing groups in fmla or letting items of `groups` be little forumlae

})

