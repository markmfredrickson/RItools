################################################################################
# Testing the behind the scenes engine 
################################################################################

library(testthat)

context("Test statistics, functions and objects")

test_that("Harmonic Mean Difference", {
  
  set.seed(20111018)
  n <- 200 # pick a bigger n than any of our other tests, as this should be fast
  Z <- rep(c(0,1), n/2)
  B <- rep(1:4, each =  n/4)
  ys <- rnorm(n) + Z + B/8 # small block effect, larger treatment effect
  
  hmd <- xBalance(Z ~ ys, data = data.frame(Z, ys), strata = as.factor(B), report =
    "adj.mean.diffs")$results[,1,1]

  expect_equal(harmonic.mean.difference(ys, Z, B), hmd)

  # the model, a location shift (does not fit reality as it ignores block
  # effects, but that's fine. no model is perfect.)
  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1)}
  
  # samples should be ignored, set low to keep the test short if there is an error
  res.xb <- randomizationDistributionEngine(ys, Z, 
                                            list(xb = list(harmonic.mean.difference, tau1)), 
                                            blocks = B,
                                            samples = 1,
                                            summaries = "z.scores",
                                            type = "asymptotic") 
  dst <- res.xb$xb

  expect_equal(dim(dst), c(2,3))
  expect_equal(colnames(dst), c("statistic", "p.value", "z"))

  # now use the exact version and check our results
  # by default type = "exact"
  res.exact <- randomizationDistributionEngine(ys, Z, 
                                            list(xb = list(harmonic.mean.difference, tau1)), 
                                            blocks = B,
                                            samples = 1)

  dst.x <- res.exact$xb

  expect_equal(dst$statistic, dst.x$statistic)

})

test_that("Wilcox.test", {
  
  set.seed(20111018)
  n <- 200 # pick a bigger n than any of our other tests, as this should be fast
  Z <- rep(c(0,1), n/2)
  B <- rep(1:4, each =  n/4)
  ys <- rnorm(n) + Z + B/8 # small block effect, larger treatment effecto
  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1) }
  
  res.wt <- randomizationDistributionEngine(ys, 
                                            Z, 
                                            list(wt = list(mann.whitney.u, tau1)),  
                                            samples = 1,
                                            type = "asymptotic")

  dst <- res.wt$wt

  expect_equal(dim(dst), c(2,2))
  expect_equal(colnames(dst), c("statistic", "p.value"))

  # the exact version
  res.exact <- randomizationDistributionEngine(ys, 
                                            Z, 
                                            list(wt = list(mann.whitney.u, tau1)),  
                                            samples = 1)

  dst.x <- res.exact$wt
  expect_equal(dst$statistic, dst.x$statistic)

})
