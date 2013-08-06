################################################################################
# Tests for xBalance function
################################################################################

library("testthat")

context("xBalance Functions")

test_that("xBalance returns covariance of tests", {
  set.seed(20130801)
  n <- 500

  library(MASS)
  xs <- mvrnorm(n, 
                mu = c(1,2,3),
                Sigma = matrix(c(1, 0.5, 0.2, 
                                 0.5, 1, 0,
                                 0.2, 0, 1), nrow = 3, byrow = T))

  p <- plogis(xs[,1]- 0.25 * xs[,2] - 1)
  z <- rbinom(n, p = p, size = 1)
  s <- rep(c(0,1), each = n/2)

  dat <- cbind(z, xs, s)


  res <- xBalance(z ~ . - s, 
                  data = as.data.frame(dat), 
                  report = 'all',
                  strata = list("Unadj" = NULL,
                                "Adj"   = ~ s))

  tcov <- attr(res$overall, "tcov")

  expect_false(is.null(tcov))

  expect_equal(length(tcov), 2)
  expect_equal(dim(tcov[[1]]), c(3,3))

  # variance should be the squares of the reported null SDs
  expect_equal(sqrt(diag(tcov[[1]])), res$results[, "adj.diff.null.sd", 1]) 
  expect_equal(sqrt(diag(tcov[[2]])), res$results[, "adj.diff.null.sd", 2]) 
})

