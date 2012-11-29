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

  # now, re-run xbalance using only X and Y to see if the "All" results match the "XY" results from previous
  xb.xy <- xBalance(Z ~ X + Y, data = df, strata = data.frame(factor("none"), df$S))
  expect_identical(xb.grp$groups[,,"XY"], xb.xy$groups[,,"All"])

  # using regular expressions to capture larger numbers of variables
  xb.any_x <- xBalance(Z ~ X + exp(X) + Y, data = df, groups = list("Xes" = c("X")))
  xb.x_expx <- xBalance(Z ~ X + exp(X), data = df)
  expect_identical(xb.any_x$groups[,,"Xes"], xb.x_expx$groups[,,'All'])

  # using regular expressions to capture categorical contributions
  # also, just use the variable name to name the gorup
  xb.all_k <- xBalance(Z ~ K, data = df, groups = list("K"))
  expect_identical(xb.all_k$groups[,,"K"], xb.all_k$groups[,,"All"])

  # double check that a group with multiple names will be pasted together
  xb.pasted_groups <- xBalance(Z ~ X * Y, data = df, groups = list(c("X", "Y")))
  expect_identical(dimnames(xb.pasted_groups$groups)[[3]], c("All", "X, Y"))

  # excluding log(X)
  xb.xonly <- xBalance(Z ~ X, data = df)
  xb.no_expx <- xBalance(Z ~ X + exp(X), data = df, groups = list("No log(X)" = "^X$"))
  expect_identical(xb.no_expx$groups[,,"No log(X)"], xb.xonly$groups[,,"All"])

  # Error checking
  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("XW" = c("W"))),
               "did not match any variables")

  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("XZ" = c("X", "Z"))),
               "Treatment")

  expect_error(xBalance(Z ~ X + Y, data = df, groups = list("Empty" = c())),
               "Empty group")


})

test_that("Creating model matrices", {
  # create a quick xBalance result
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = rep(c(T,F), each = n/2),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   S = as.factor(rep(c("A", "B"), each = n/2)))

  f <- terms(Z ~ X * K + W + S, data = df)

  mm <- make_nice_model_matrix(f, model.frame(f, df))

  expect_is(mm, "matrix")
  expect_true(!("(Intercept)" %in% colnames(mm)))
  expect_equal(dim(mm), c(100, 9)) # 3 for each K level, 3 for X * K interactions, 1 for X, 1 for W = TRUE, 1 for S = B
  expect_true(all(c("K = a", "K = b", "K = c", "X:K = a", "X:K = b", "X:K = c", 
                    "X", "W", "S") %in% colnames(mm)))
  
  # the assignnames attribute should match the columns to the name of each variable in the original formula/model.frame
  # this will be used to create and use groups
  expect_equal(length(unique(attr(mm, "assignnames"))), 5)
  vars <- c("X", "K", "W", "S", "X:K")
  assignnames <- attr(mm, "assignnames")
  expect_true(all(vars %in% assignnames) && all(assignnames %in% vars))


  # make sure that results are matrices, even when using a single var
  mm.1 <- make_nice_model_matrix(Z ~ X, model.frame(Z ~ X, df))
  expect_is(mm.1, "matrix")
  expect_equal(dim(mm.1), c(100,1))

  # NAs should be allowed to pass through
  df.NA <- df
  df.NA[1:5, "X"] <- NA
  
  mf.NA <- model.frame(Z ~ X, data = df.NA, na.action = na.pass)
  expect_equal(dim(mf.NA), c(100, 2))

  mm.NA <- make_nice_model_matrix(terms(Z ~ X), mf.NA)
  expect_equal(dim(mm.NA), c(100, 1))

  # expressions in the LHS and RHS
  f.lhs <- terms(I(W > 0) ~ X)
  f.rhs <- terms(Z ~ I(X^2))

  mf.lhs <- model.frame(f.lhs, df)
  mf.rhs <- model.frame(f.rhs, df)

  mm.lhs <- make_nice_model_matrix(f.lhs, mf.lhs)
  mm.rhs <- make_nice_model_matrix(f.rhs, mf.rhs)

  expect_equal(colnames(mm.lhs), "X")
  expect_equal(colnames(mm.rhs), "I(X^2)")

})

