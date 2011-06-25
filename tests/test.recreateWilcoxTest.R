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
  # the result of the test stat applied to the observed data it the first
  # entry in the sharp null object
  expect_equal(res.prd[[1,1]], res.wilcox$statistic[[1]])

  s <- summary(res.prd)
  expect_equal(res.wilcox$p.value, s@sharp.null.p)
  
  res.wilcox.m10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = -10)
  res.wilcox.10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 10)
  res.wilcox.9 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 9)

  pvs <- na.omit(p.values(res.prd)) # na.omit drops the sharp null, which complicates
  
  expect_equal(pvs[pvs$tau == 10, "p"], res.wilcox.10$p.value)
  expect_equal(pvs[pvs$tau == -10, "p"], res.wilcox.m10$p.value)
  expect_equal(pvs[pvs$tau == 9, "p"], res.wilcox.9$p.value)

  # zeroing in on the point estimate
  step <- 0.01
  res.prd.intval <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(9,10.5, step)))

  # this problem does not have a unique point estimate. Typical HL would then take the
  # mid point of the interval. In our case, we'll take the mean and call it close enough.
  s.intval <- summary(res.prd.intval)
  thept <- mean(s.intval@point.estimate$tau)
  expect_true(abs(thept - res.wilcox$estimate[[1]]) < step) # [[]] to remove name
})

