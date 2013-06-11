################################################################################
# Tests comparing models over various distance metrics
################################################################################

library("testthat")

context("Model Distance Comparisons")

test_that("parameterSensitivity", {
  # a function that takes a model, a set of params, and compares the distance
  # of prediction versus the distance of parameters

  expect_true(exists("parameterSensitivity"))

  # the result should be a matrix with two columns (prediction distance, param
  # distance) with one row for every pairwise comparison of models

  # start with a dummy model that returns the same value for all units, always
  # but has two meaningless parameters, a and b
  const.model <- UniformityModel(
    function(y, z, a, b) { y }, 
    function(y, z, a, b) { y })

  res <- parameterSensitivity(const.model,
                              parameters = list(a = 1:5, b = -3:7),
                              uniformity = runif(100),
                              z = rep(c(0,1), 50))

  na <- 5
  nb <- 11
  nab <- 5 * 11

  expect_equal(dim(res), c(nab * (nab - 1) / 2, 6))
  expect_equal(colnames(res), c("left.a", "left.b", "right.a", "right.b", "parameter", "prediction"))
  expect_true(all(is.nan(res[,"prediction"])))

  # now try with a slightly more interesting model
  res.c <- parameterSensitivity(constant.additive.model,
                                parameters = list(tau = 0:20),
                                uniformity = rnorm(100, 10),
                                z = rep(c(0,1), 50))

  expect_true(!all(res.c[,"prediction"] == 0))
  expect_equal(colnames(res.c), c("left.tau", "right.tau", "parameter", "prediction"))
  expect_true(all(res.c[, c("parameter", "prediction")] <= 1))
  expect_true(all(res.c[, c("parameter", "prediction")] >= 0))
  
  # select an intentionally insensitive model
  double.add <- UniformityModel(
    function(y, z, a, b) { y - z * a - z * b },
    function(y, z, a, b) { y + z * a + z * b })

  res.double.add <- parameterSensitivity(double.add,
                                         parameters = list(a = 1:10, b = -5:5),
                                         uniformity = rnorm(100, 10),
                                         z = rep(c(0,1), 50))

  expect_true(all(res.double.add[, c("parameter", "prediction")] <= 1))
  expect_true(all(res.double.add[, c("parameter", "prediction")] >= 0))

  
})

