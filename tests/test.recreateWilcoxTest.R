################################################################################
# Testing the parameterizedRandomizationDistribution against wilcox.test
################################################################################

library(testthat)

context("Wilcox-like results")

test_that("Get same point estimate", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  expect_equal(sum(Z), 6)

  res.wilcox <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T)

  res.prd <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = c(-10, 9, 10)))

  # the wilcox stat is named "W", but is otherwise the same
  expect_equal(res.prd@observed.test.stat, res.wilcox$statistic[[1]])

  s <- summary(res.prd)
  expect_equal(res.wilcox$p.value, s@sharp.null.p)
  
  res.wilcox.m10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = -10)
  res.wilcox.10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 10)
  res.wilcox.9 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 9)

  pvs <- p.values(res.prd)
  
  expect_equal(pvs[pvs$tau == 10, "p"], res.wilcox.10$p.value)
  expect_equal(pvs[pvs$tau == -10, "p"], res.wilcox.m10$p.value)
  expect_equal(pvs[pvs$tau == 9, "p"], res.wilcox.9$p.value)

})

