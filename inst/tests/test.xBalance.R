################################################################################
# Tests for xBalance function
################################################################################
## if working interactively in inst/tests you'll need
## library(RItools, lib.loc = '../../.local')
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


  # we use ETT weighting here to correspond to the weighting scheme used
  # in the descriptives section
  res <- xBalance(z ~ . + strata(s),
                  data = as.data.frame(dat),
                  stratum.weights = RItools:::effectOfTreatmentOnTreated,
                  report = 'all')

  tcov <- attr(res$overall, "tcov")

  expect_false(is.null(tcov))

  expect_equal(length(tcov), 2)
  expect_equal(dim(tcov[[1]]), c(4,4))

})

test_that("Passing post.alignment.transform, #26", {
  data(nuclearplants)

  # Identity shouldn't have an effect
  res1 <- xBalance(pr ~ ., data=nuclearplants)
  res2 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = function(x) x)

  expect_true(all.equal(res1, res2)) ## allow for small numerical differences

  res3 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = rank)

  expect_true(all(dim(res1$results) == dim(res3$results)))

  expect_error(xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = mean),
               "Invalid post.alignment.transform given")

  res4 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = rank, report="all")
  res5 <- xBalance(pr ~ ., data=nuclearplants, report="all")

  expect_false(isTRUE(all.equal(res4,res5)))

  # a wilcoxon rank sum test, asymptotic and w/o continuity correction
  res6 <- xBalance(pr ~ cost, data=nuclearplants, post.alignment.transform = rank, report="all")

  expect_equal(res6$results["cost", "p", "Unstrat"],
               wilcox.test(cost~pr, data=nuclearplants, exact=FALSE, correct=FALSE)$p.value)

  # w/ one variable, chisquare p value should be same as p value on that variable
  expect_equal(res6$results["cost", "p", "Unstrat"],
               res6$overall["Unstrat","p.value"])

  # to dos: test combo of a transform with non-default stratum weights.

})

test_that("NA in stratify factor are dropped", {
  data(nuclearplants)

  n2 <- nuclearplants
  n2 <- rbind(n2, n2[1,])
  n2$pt[1] <- NA

  f <- function(d) {
    xBalance(pr ~ . - pt + strata(pt) - 1, data = d)
  }

  xb1 <- f(nuclearplants)
  xb2 <- f(n2)

  expect_equal(xb1, xb2)
})
