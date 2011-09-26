################################################################################
# Objects and functions to work with models
################################################################################

library(testthat)
context("Randomization Model Tests")

test_that("Basics", {
  # tau helper returns the passed in tau value:
  expect_equal(returnTau(beta = 1, 3, 4, tau = 5, 7, gamma = -3), 5)

  # test the constant additive model
  uniformity <- c(0, 1, 2, 4)
  Z <- c(0, 1, 0, 1)
  true.tau <- 3
  expected <- c(0, 4, 2, 7)
  
  expect_equal(
    modelOfEffect(constant.additive.model, expected, Z, blocks = NULL, tau =
      true.tau),
    uniformity)
  
})

test_that("Network models", {
  # super duper simple network model
  # get an effect of 1 if any neighbors, zero otherwise
  # network: 1 <-> 2, 2 <-> 3, 4 has no neighbors
  S <- matrix(c(0,1,0,0, 1,0,1,0, 0,1,0,0, 0,0,0,0), nrow = 4)
  uniformity <- c(0,1,2,4)
  Z <- c(0, 1, 0, 0)
  true.tau <- 3
  expected.output <- c(1, 4, 3, 4)
  sdsn.model <- networkRandomizationModel(S, returnTau, 
    function(Z, blocks, ZS, tau) {
      ifelse(ZS > 0, 1, 0)  
    })

  expect_equal(modelOfEffect(sdsn.model, expected.output, Z, blocks = NULL, tau = true.tau),
    uniformity)

  expect_equal(observedData(sdsn.model, uniformity, Z, blocks = NULL, tau = true.tau),
    expected.output)

})

test_that("functions", {
  # functions should work as modelsOfEffect and observedData first arguments
  constant.additive.moe <- function(ys, z, b, tau) {
    ys - (z * tau)
  } 
  constant.additive.obdata <- function(ys, z, b, tau) {
    ys + (z * tau)
  } 
  
  data <- 1:6
  Z <- rep(c(1,0), 3)
  true.tau <- 4

  expect_equal(
    modelOfEffect(constant.additive.model, data, Z, NULL, tau = true.tau),
    modelOfEffect(constant.additive.moe, data, Z, null, tau = true.tau))

  expect_equal(
    observedData(constant.additive.model, data, Z, NULL, tau = true.tau),
    observedData(constant.additive.obdata, data, Z, null, tau = true.tau))

  
})

test_that("Model Analysis Functions", {
  
  Z <- c(1,0,1,0)
  data <- c(1,2,3,4)

  res.analysis <- analyzeModel(constant.additive.model,
                               Z,
                               list(tau = 3), 
                               list(tau = 0:5),
                               mean.diff.noblocks)

  expect_equal(length(res.analysis$simulation), 6)

  

})

test_that("Error checking", {
  Z <- c(1,0,1,0)
  data <- c(1,2,3,4)

  expect_error(analyzeModel(constant.additive.model,
                               Z,
                               list(tau = 3), 
                               list(tau = 0:5, beta = 1),
                               mean.diff.noblocks))
  
})
