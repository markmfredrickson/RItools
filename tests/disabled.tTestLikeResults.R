################################################################################
# Testing the RItest against t.test
################################################################################

library(testthat)

context("t Test like results")

test_that("Get same point estimate", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  res.t.test <- t.test(R[Z == 1], R[Z == 0], var.equal = T)
  
  res.prd <- RItest(R, Z, mean.difference, 
    constant.additive.model, list(tau = c(-10, 9, 10)))
  
  # the t.test p-value is obscenely small, so just make sure that the null
  # pvalue is equal to the smallest real possible p (two tailed)
  minp <- 2 * 1/choose(n, n/2)
  expect_equal(res.prd[1,2], minp)

  # this tolerance value is arbitrary, but there you are.
  within <- function(a, b, tol = 0.01) {
    expect_true(abs(a - b) <= tol)  
  }

  res.t.test.10 <- t.test(R[Z == 1], R[Z == 0], var.equal = T, mu = 10)
  res.t.test.9 <- t.test(R[Z == 1], R[Z == 0], var.equal = T, mu = 9)

  within(res.prd[4,2], res.t.test.10$p.value)
  within(res.prd[3,2], res.t.test.9$p.value)
})
