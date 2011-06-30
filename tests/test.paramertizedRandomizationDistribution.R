################################################################################
# Test of the pRD() function. Many other tests of this function in other test.*
# files. This is a place for tests that don't have a home there.
################################################################################

library(testthat)
context("paramertizedRandomizationDistribution")

test_that("Engine and pRD give same answers", {
  
  # data borrowed from wilcox test
  set.seed(20110620)
  tau <- 10
  n <- 8 
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  set.seed(20110620)
  res.prd <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = c(-10, 9, 10)))

  # being extra explicit about model creation
  mm10 <- function(y, z, b) { constant.additive.model(y, z, b, -10) }
  m9 <- function(y, z, b) { constant.additive.model(y, z, b, 9) }
  m10 <- function(y, z, b) { constant.additive.model(y, z, b, 10) }
  
  set.seed(20110620)
  res.eng <- randomizationDistributionEngine(R, Z, list(list(mann.whitney.u, mm10, m9, m10)))

  # prd has some extra information, and we don't expect that to be the same
  # also, the return value of res.eng is a list, so we pull out the first item, just as pRD does.
  expect_identical(res.prd@.Data, res.eng[[1]]@.Data)
})

