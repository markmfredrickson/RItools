################################################################################
# Tests of randomziation producing functions
################################################################################

library(testthat)

context("Randomization functions")

test_that("Correct number of Zs", {

  ### Check bad input
  # Z and B must be same length
  expect_error(produceRandomizations(rep(1, 10), rep(1, 5), 10))

  ### Proper Sample Sizes
  ## Small problems
  # choose(10, 5) => 252
  res.noblks <- produceRandomizations(rep(c(0,1), 5), rep(1, 10), samples =
    1000)
  expect_equal(dim(res.noblks), c(5, 252))

  res.noblks2 <- produceRandomizations(rep(c(0,1), 5), rep(1, 10), samples =
    251)
  expect_equal(dim(res.noblks2), c(5, 251))

  res.noblks3 <- produceRandomizations(rep(c(0,1), 5), rep(1, 10), samples =
    253)
  expect_equal(dim(res.noblks3), c(5, 252))

  # 2^5 = 32
  res.blkd <- produceRandomizations(rep(c(0,1), 5), rep(1:5, 2), samples = 1000)
  expect_equal(dim(res.blkd), c(5, 32))

  ## Big problems
  res.big <- produceRandomizations(rep(c(0,1), 100), rep(1, 200), samples =
    500)
  expect_equal(dim(res.big), c(100, 500))

})
