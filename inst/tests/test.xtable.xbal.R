################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("xtable methods")

test_that("xtable xbal method", {
  set.seed(20121129)

  # create a quick xBalance result
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = cut(rnorm(n), breaks = 3),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   U = rnorm(n),
                   S = as.factor(rep(c("A", "B"), each = n/2)))


  # the data are will be grouped as (X, XY) and (all levels W by all levels
  # K), Y is not in an explicit group
  grps <- list("X (with interactions)" = "X", 
               "Categorical Variables" = c("W", "K"),
               "Interactions" = ":")

  xb <- xBalance(Z ~ X * Y + W * K, 
                 data = df, 
                 strata = data.frame(none = factor("none"), other = df$S),
                 groups = grps)

  # just double check some properties of the xbal object
  # df = X:Y + (W:K => 3 * 3 - 1)
  expect_equal(xb$groups[1, "df", "Interactions"], 9)

  library(xtable)
  xt <- xtable(xb)

  expect_is(xt, "xtable")

  # an xtable object is basically a data.frame
  # we should expect one row for each variable 
  # X, Y, X:Y, 3 * W, 3 * K, 3 * 3 (W * K)
  expect_equal(dim(xt)[1], 18)

  xt.grp <- xtable(xb, type = "groups")
  expect_true(!identical(xt, xt.grp))

  expect_is(xt.grp, "xtable")

  # like the variables table, it is a flattened version of an array, with rows
  # of groups, columns of strata by tests
  expect_equal(dim(xt.grp), c(4, 6))
  expect_equal(rownames(xt.grp), c("All", names(grps)))
  expect_equal(colnames(xt.grp), rep(c("chisquare", "df", "p.value"), 2))


  # check that the default is just the variables
  expect_identical(xt, xtable(xb, type = "variables"))

  # error checking
  expect_error(xtable(xb, type = "foo"), "Unknown type of table: foo")
})

