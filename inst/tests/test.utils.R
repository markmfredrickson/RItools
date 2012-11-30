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

  # we expect to have a no intercept formula
  expect_equal(Z ~ X, as.formula(xb))
})

test_that("Select a subset of xbal results (for printing, etc)", {
  # create a quick xBalance result
  set.seed(20121129)
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
  xb <- xBalance(Z ~ X * Y + W * K, 
                 data = df, 
                 strata = data.frame(none = factor("none"), 
                                     S = df$S, 
                                     U = cut(df$U, 3)),
                 groups = list("X (with interactions)" = "X", 
                               "Categorical Variables" = c("W", "K"),
                               "Interactions" = ":"))

  # just double check some properties of the xbal object
  # df = X:Y + (W:K => 3 * 3 - 1)
  expect_equal(xb$groups[1, "df", "Interactions"], 9)

  # strata based selection is the easiest as it common across groups and
  # variables
  xb.noneU <- select(xb, which.strata = c("none", "U"))
  
  expect_equal(dim(xb$results)[["strata"]], 3)
  expect_equal(dim(xb$groups)[["strata"]], 3)
  expect_equal(dim(xb.noneU$results)[["strata"]], 2)
  expect_equal(dim(xb.noneU$groups)[["strata"]], 2)
   

})

