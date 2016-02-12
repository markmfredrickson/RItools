library(testthat)
context("RItest")

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
  
  expect_is(res.model, "array")
  expect_equal(dim(res.model), sapply(hypotheses, length))

  model3 <- function(y, z, a, b, c) {
    y - z * a * (b ^ c)
  }

  hypo3 <- list(a = 1:10, b = -5:5, c = c(2,4,6))

  res3 <- RItest(R, Z, mean.difference, model3, parameters = hypo3, samples = 10)

  expect_equal(dim(res3), sapply(hypo3, length))
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

  model2 <- function(y, z, tau, gamma) {
    ifelse(!!z, y * gamma - tau, y)  
  }

  res.two <- RItest(R, Z, mean.difference,
    model2, parameters = hypotheses)

  a <- plot(res.two)
  expect_equal(a$panel, "panel.levelplot")

})
