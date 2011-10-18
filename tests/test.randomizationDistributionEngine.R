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

  res <- randomizationDistributionEngine(ys, Z, models)

  # res is a list of results, same length as models
  expect_is(res, "list")
  expect_equal(length(res),1)
  expect_equal(names(res), c("xyz"))

  # each list time should be a RandomizationDistribution object
  dst <- res$xyz
  expect_is(dst, "RandomizationDistribution")
  
  expect_true(all(sapply(c("test.statistic", "models.of.effect", "treatment", "blocks"),
                         function(s) { !(is.null(slot(dst, s)))})))

  # dst inherits from a matrix, and should have dim = c(models, samples + 1)
  # two models were tested (sharp null and tau1), there are choose(8,4) = 70 samples
  # the extra in the sample column is for the observed test statistics (after adjustment)
  # this case, there should be another value in the table with exactly the same value as we are enumerating.
  expect_equal(dim(dst), c(2, 71))
  
  
})
