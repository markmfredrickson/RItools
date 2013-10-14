################################################################################
# Testing the RItest against wilcox.test
################################################################################

library(testthat)

context("Wilcox-like results")

test_that("Non-paired", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  expect_equal(sum(Z), 6)

  res.wilcox <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, alternative = "greater")

  res.prd <- RItest(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = c(-10, 9, 10)))

  # the wilcox stat is named "W", but is otherwise the same
  # the result of the test stat applied to the observed data it the first
  # entry in the sharp null object
  expect_equal(res.prd@observed.statistic, res.wilcox$statistic[[1]])

  expect_equal(res.wilcox$p.value, as.numeric(res.prd@sharp.null))
  
  res.wilcox.m10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = -10, alternative = "greater")
  res.wilcox.10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 10, alternative = "greater")
  res.wilcox.9 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 9, alternative = "greater")

  expect_equivalent(as.numeric(res.prd), c(res.wilcox.m10$p.value, 
      res.wilcox.9$p.value, res.wilcox.10$p.value))
})

test_that("Paired", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  nB <- n/2 ## number of blocks, here, pairs
  
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(c(0,1), n/2) ; z <- as.logical(Z) # for convenience, create a logical version
  R <- Z * Yt + (1 - Z) * Yc
  B <- gl(nB,2)
  
  expect_equal(sum(Z), 6)
  res.wilcox <- wilcox.test(R[z], R[!z], paired = T, exact = T, conf.int = T, alternative = "greater")

  # make sure the test statistic is correct on the observed data
  expect_equal(paired.sgnrank.sum(R, Z), res.wilcox$stat[[1]])

  res.prd <- RItest(R, Z, paired.sgnrank.sum, 
    constant.additive.model, list(tau = c(-10, 9, 10)),
    sampler = simpleRandomSampler(z = Z, b = B)) 

  res.wilcox.m10 <- wilcox.test(R[z], R[!z], paired = T, exact = T, conf.int = T, mu = -10, alternative = "greater")
  res.wilcox.10 <-  wilcox.test(R[z], R[!z], paired = T, exact = T, conf.int = T, mu = 10, alternative = "greater")
  res.wilcox.9 <-   wilcox.test(R[z], R[!z], paired = T, exact = T, conf.int = T, mu = 9, alternative = "greater")

  expect_equivalent(as.numeric(res.prd), c(res.wilcox.m10$p.value, 
      res.wilcox.9$p.value, res.wilcox.10$p.value))

})

