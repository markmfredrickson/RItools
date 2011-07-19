################################################################################
# Attributable effects tests
################################################################################  

library(testthat)

context("Attributable Effects")

test_that("Helpers", {
  n <- 12
  Yc <- c(rep(1, n/3), rep(0, 2 * n/3))
  Yt <- c(rep(1, 2 * n/3), rep(0, n/3)) 
  Z <- rep(c(0,1), n/2)
  R <- Z * Yt + (1 - Z) * Yc

  expect_equal(max.attributable.effects(R, Z), 4)
  expect_equal(min.attributable.effects(R, Z), -2)
  
  # model gets the numbers right, though not the exact units (which are
  # exchangeable anyway)
  expect_equal(sum(additive.attributable.effect(R, Z, NULL, sum(R - Yc))),
    sum(Yc))

  expect_equal(sum(additive.attributable.effect(R, Z, NULL, 1)) + 1, sum(R))
  expect_equal(sum(additive.attributable.effect(R, Z, NULL, -1)) - 1, sum(R))

})

test_that("Checks input", {
  # outcome must be binary: 0 or 1
  Y_bad1 <- rep(c(1,2), 5)
  Y_bad2 <- rep(0:4, 2)
  Z <- rep(c(0,1), 5)

  expect_error(attributableEffects(Y_bad1, Z, mean.difference))
  expect_error(attributableEffects(Y_bad2, Z, mean.difference))

})
