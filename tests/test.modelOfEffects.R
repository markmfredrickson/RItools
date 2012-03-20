################################################################################
# Tests for models of effects and related functions
################################################################################

library(testthat)

context("Small models and helper functions")

test_that("givenParams helper", {
  mymodel <- function(y, z, b, q, h) { y - z * q + h }
  Z <- c(1,0,0,1)
  Y <- c(1,2,3,4)
  
  mymodel.two <- givenParams(mymodel, q = 2)
  expect_equal(mymodel.two(Y, Z, NULL, h = 1), c(0,3,4,3))

})

test_that("min.max.model meta model", {
  additive <- function(y, z, b, tau) { y - z * tau}
  Z <- c(1,0,0,1)
  Y <- c(1,2,3,4)
  
  additive.0.3 <- min.max.model(additive, 0, 3)

  expect_equal(additive.0.3(Y, Z, NULL, tau = 3), c(0,2,3,1))
})

