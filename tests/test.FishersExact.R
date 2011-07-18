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
  
  lady.distribution <- parameterizedRandomizationDistribution(lady.guess, actual.cups,
    test.stat = number.correct)

  pv <- upper.p.value(lady.distribution[1,1],
    lady.distribution[1,-1])
  expect_equal(pv, 17/70)

  expect_equivalent(as.numeric(table(lady.distribution[1,-1])/70), 
    c(1/70, 16/70, 36/70, 16/70, 1/70))
  
})

test_that("2x2 table style", {
  set.seed(20110620)
  n <- 12
  Yc <- c(rep(1, n/3), rep(0, 2 * n/3))
  Yt <- c(rep(1, 2 * n/3), rep(0, n/3)) # added 4 units
  Z <- rep(c(0,1), n/2)
  R <- Z * Yt + (1 - Z) * Yc
  
  res <- parameterizedRandomizationDistribution(R, Z, test.stat = odds.ratio)

  res.fisher <- fisher.test(table(R, Z))

  res.p <- general.two.sided.p.value(res[1,1], res[1,-1])
  
  expect_equal(res.p, res.fisher$p.value)


})

