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

  # dst inherits from a matrix, and should have dim = c(models, 2)
  # two models were tested (sharp null and tau1), there are choose(8,4) = 70 samples
  # the first column is the observed statistic after adjustment, the second the p-value
  # computed using the p.value argument (usually general.two.sided.p.value)
  expect_equal(dim(dst), c(2, 2))
  expect_equal(dst@samples, 70)

  # since we didn't supply a p.value function, expect it equal to general.two.sided.p.value
  expect_equal(general.two.sided.p.value, dst@p.value)

  # blocks should just be a vector of 1s, as we didn't supply any blocks
  expect_equal(dst@blocks, rep(1, 8))
  
})

test_that("Get entire distribution", {
  Z <- rep(c(0,1), 4)
  ys <- rnorm(8) + Z
  
  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1)}

  models <- list(xyz = list(mean.difference, tau1)) # tests sharp null and tau1

  # typically we don't return the entire distribution, only the p-values
  res <- randomizationDistributionEngine(ys, Z, models, include.distribution = TRUE)
 
  dst <- res$xyz
  
  # 2 models, choose(8, 4) = 70 randomizations
  expect_equal(dim(dst@distribution), c(2, 70))
  
})

test_that("Summaries of the distrib", {
  Z <- rep(c(0,1), 4)
  ys <- rnorm(8) + Z
  
  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1)}

  models <- list(xyz = list(mean.difference, tau1)) # tests sharp null and tau1

  res <- randomizationDistributionEngine(ys, Z, models,
            summaries = list(mean = mean, var = var))

  dst <- res$xyz
  expect_equal(dim(dst), c(2,4))
  expect_equal(colnames(dst), c("statistic", "p.value", "mean", "var"))
    
})

test_that("Multiple backends", {
  # the engine has two steps: 
  # - a front end that adjusts data consistent with a model
  # - a backend that computes the null distribution of this data
  # typically we use the re-sampling backend that shuffles Z
  # we could use xBalance and Wilcox test as well, as they are valid
  # randomization tests, provided we provide them the adjusted data
  # this section includes tests for this functionality.


  # xBalance and wilox.test are a test statistic of central tendancy, but
  # that does not imply the model need be a location shift. 
  # we will use a location shift, though it could be any model

  n <- 200 # pick a bigger n than any of our other tests, as this should be fast
  Z <- rep(c(0,1), n/2)
  B <- rep(1:4, each =  n/4)
  ys <- rnorm(n) + Z + B/8 # small block effect, larger treatment effect

  # the model, a location shift
  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1)}

  # first, let's test if we can get a backend called.
  # any object with the "randomizationEngine" attr calls the attribute
  # this is how we'll extend exising functions (e.g. xBalance and wilcox.test)
  
  test.backend <- function(...) { stop("Backend called") }
  test.backend.test.stat <- "hello!"
  attr(test.backend.test.stat, "randomizationEngine") <- test.backend

  expect_error(randomizationDistributionEngine(ys, Z, list(xb = list(test.backend.test.stat, tau1))),
    "Backend called")

  # now on to the fun stuff!
  # samples should be ignored, set low to keep the test short if there is an error
  res.xb <- randomizationDistributionEngine(ys, Z, list(xb = list(xBalance, tau1)), 
                                            blocks = B,
                                            samples = 1,
                                            summaries = "z.scores") 
 
  dst <- res.xb$xb

  expect_equal(dim(dst), c(2,3))
  expect_equal(colnames(dst), c("statistic", "p.value", "z"))

  # now for the wilcox.test backend
  res.wt <- randomizationDistributionEngine(ys, Z, list(wt = list(wilcox.test, tau1)),  
                                            samples = 1)

  dst <- res.wt$wt

  expect_equal(dim(dst), c(2,2))
  expect_equal(colnames(dst), c("statistic", "p.value"))

  

})
