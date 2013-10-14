################################################################################
# Fisher's Exact test, with classic examples
################################################################################

library(testthat)
context("Fisher exact test")

test_that("Lady Tasting Tea", {
  actual.cups <- c(1, 0, 0, 1, 1, 0, 1, 0)
  lady.guess <- c(1, 0, 0, 1, 1, 0, 0, 1)

  number.correct <- function(guesses, z, blocks) { sum(z == guesses) / 2 }
  expect_equal(number.correct(lady.guess, actual.cups, NULL), 3)
  
  lady <- RItest(lady.guess, actual.cups,
    test.stat = number.correct)

  expect_equal(as.numeric(lady), 17/70)

})

test_that("2x2 table style", {
  set.seed(20110620)
  n <- 12
  Yc <- c(rep(1, n/3), rep(0, 2 * n/3))
  Yt <- c(rep(1, 2 * n/3), rep(0, n/3)) # added 4 units
  Z <- rep(c(0,1), n/2)
  R <- Z * Yt + (1 - Z) * Yc
  
  res <- RItest(R, Z, test.stat = odds.ratio)

  res.fisher <- fisher.test(table(R, Z), alternative = "greater")
  
  expect_equal(as.numeric(res), res.fisher$p.value)


})

