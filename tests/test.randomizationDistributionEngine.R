################################################################################
# Testing the behind the scenes engine 
################################################################################

library(testthat)

context("Randomization Distribution Engine")

test_that("Basic usage", {
  Z <- rep(c(0,1), 4)
  ys <- rnorm(8) + Z
  
  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1)}

  models <- list(xyz = list(mean.difference, tau1)) # tests sharp null and tau1

  res <- randomizationDistributionEngine(ys, Z, models, samples = 10000) # there will only be 70 samples

  # res is a list of results, same length as models
  expect_is(res, "list")
  expect_equal(length(res),1)
  expect_equal(names(res), c("xyz"))

  # each list time should be a RandomizationDistribution object
  dst <- res$xyz
  expect_is(dst, "RandomizationDistribution")
  
  # this test is pretty minimal as R automatically fills in non-null values (such as empty vectors)
  expect_true(all(sapply(c("test.statistic", "models.of.effect", "treatment", "blocks", "samples", "p.value"),
                         function(s) { !(is.null(slot(dst, s)))})))

  # these should all have at least one value, more for treatment and blocks
  expect_true(all(sapply(c("treatment", "blocks", "samples", "models.of.effect"),
                         function(s) { length(slot(dst, s)) > 0 })))

  # the "null" function supplied by R automatically has no arguments. These should all take a few args.
  expect_true(all(sapply(c("test.statistic", "p.value"),
                         function(s) { length(formals(slot(dst, s))) > 0 })))

  # dst inherits from a matrix, and should have dim = c(models, samples + 1)
  # two models were tested (sharp null and tau1), there are choose(8,4) = 70 samples
  # the extra in the sample column is for the observed test statistics (after adjustment)
  # this case, there should be another value in the table with exactly the same value as we are enumerating.
  expect_equal(dim(dst), c(2, 71))
  expect_equal(dst@samples, 70)

  # since we didn't supply a p.value function, expect it equal to general.two.sided.p.value
  expect_equal(general.two.sided.p.value, dst@p.value)

  # blocks should just be a vector of 1s, as we didn't supply any blocks
  expect_equal(dst@blocks, rep(1, 8))
  
  
})
