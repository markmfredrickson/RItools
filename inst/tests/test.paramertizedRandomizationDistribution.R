################################################################################
# Test of the pRD() function. Many other tests of this function in other test.*
# files. This is a place for tests that don't have a home there.
################################################################################

library(testthat)
context("paramertizedRandomizationDistribution")

test_that("Using model objects", {
  # constant.additive.model is already defined
  set.seed(20110620)
  tau <- 10
  n <- 8 
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  function.version <- function(y,z,tau) {
    y - z * tau
  }

  hypotheses <-  list(tau = c(7,8,9,10,11,12))
  res.fn <- RItest(R, Z, mean.difference, 
    function.version, parameters = hypotheses)

  res.model <- RItest(R, Z, mean.difference,
    constant.additive.model, parameters = hypotheses)
  
})

test_that("Moe can be NULL", {
  set.seed(20110620)
  tau <- 10
  n <- 8 
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  # this should not raise an error
  res.prd <- RItest(R, Z, mean.difference)
  
})

test_that("Plotting", {
  # these tests mostly just test for existence, to make sure things don't break  
  # for correct presentation, these should be tested by hand.
  
  # one param model
  set.seed(20110620)
  tau <- 10
  n <- 8 
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  hypotheses <-  list(tau = c(7,8,9,10,11,12))

  res.one <- RItest(R, Z, mean.difference,
    constant.additive.model, parameters = hypotheses)

  a <- plot(res.one)
  expect_equal(a$panel, "panel.xyplot")
  
  # two param model
   
  gamma <- 4
  Yt <- Yc/gamma + tau
  R <- Z * Yt + (1 - Z) * Yc

  hypotheses <- list(tau = 5:15, gamma = 0:5)

  model2 <- function(y, z, b, tau, gamma) {
    ifelse(!!z, y * gamma - tau, y)  
  }

  res.two <- RItest(R, Z, mean.difference,
    model2, parameters = hypotheses)

  a <- plot(res.two)
  expect_equal(a$panel, "panel.levelplot")

  # can't plot sharp nulls
  res.sn <- RItest(R, Z, mean.difference)
  expect_error(plot(res.sn), "Cannot plot sharp null only")

})
