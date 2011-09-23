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

