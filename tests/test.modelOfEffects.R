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

  # use the constant.additive.model which is a UniformityModel object
  add2 <- givenParams(constant.additive.model, tau = 2)
  expect_equal(c(0,0,0), add2(c(2,0,2), c(1,0,1)))
  expect_equal(c(2,0,2), invertModel(add2, c(0,0,0), c(1,0,1)))

})

test_that("min.max.model meta model", {
  additive <- function(y, z, b, tau) { y - z * tau}
  Z <- c(1,0,0,1)
  Y <- c(1,2,3,4)
  
  additive.0.3 <- min.max.model(additive, 0, 3)

  expect_equal(additive.0.3(Y, Z, NULL, tau = 3), c(0,2,3,1))
})

