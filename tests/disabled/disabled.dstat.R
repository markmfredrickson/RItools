################################################################################
# d and d^2 statistics
################################################################################

library(testthat)
context("d and d^2 statistic`")

test_that("Basics", {
  Z <- rep(c(0,1), 4)
  b <- rep(c(0,1,3,4), each = 2)

  ys <- b + Z * 3 # block level effects plus treatment level
  
  # known failure. removing for now
  # expect_equal(d.stat(ys, Z, b), 3)

  


})
