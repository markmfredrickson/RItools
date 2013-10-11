################################################################################
# Testing the behind the scenes engine 
################################################################################

library(testthat)

context("Randomization Distribution Engine")

test_that("Basic usage", {
  Z <- rep(c(0,1), 4)
  ys <- rnorm(8) + Z
  
  tau1 <- function(y, z) { constant.additive.model(ys, z, tau = 1)}

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
  expect_true(all(sapply(c("test.statistic", "models.of.effect", "z", "samples", "p.value"),
                         function(s) { !(is.null(slot(dst, s)))})))

  # these should all have at least one value, more for treatment 
  expect_true(all(sapply(c("z", "samples", "models.of.effect"),
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

})

test_that("Get entire distribution", {
  Z <- rep(c(0,1), 4)
  ys <- rnorm(8) + Z
  
  tau1 <- function(y, z) { constant.additive.model(ys, z, tau = 1)}

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
  
  tau1 <- function(y, z) { constant.additive.model(ys, z, tau = 1)}

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
  tau1 <- function(y, z) { constant.additive.model(ys, z, tau = 1)}
 
  # first, let's test if we can get a backend called.
  # any object with the "randomizationEngine" attr calls the attribute
  # this is how we'll extend exising functions (e.g. xBalance and wilcox.test)
  
  test.backend <- new("AsymptoticTestStatistic", 
    function(y,z) { stop("Not asymptotic") },
    asymptotic = function(...) { stop("Backend called")})

  expect_error(randomizationDistributionEngine(ys, Z, list(xb = list(test.backend, tau1)), type = "asymptotic"),
    "Backend called")

  expect_error(randomizationDistributionEngine(ys, Z, list(xb = list(test.backend, tau1)), type = "exact"),
    "Not asymptotic")
  expect_error(randomizationDistributionEngine(ys, Z, list(xb = list(test.backend, tau1))), 
    "Not asymptotic")
 
  # real errors for non asym test stats

  ts <- function(y, z) { 1 }
  expect_error(randomizationDistributionEngine(ys, Z, list(xb = list(ts)), type = "asymptotic"), "asymptotic") # use sharp null)

  # should not call the sampler
  safeBackEnd <- new("AsymptoticTestStatistic", 
    function(y,z) { return(TRUE)},
    asymptotic = function(...) { return(FALSE) }) 
  
  dummySampler <- function(n) { stop("Oops!") }
  randomizationDistributionEngine(ys, Z, list(xb = list(safeBackEnd)), type = "asymptotic", sampler = dummySampler)

  # only valid types are exact and asymptotic
  expect_error(randomizationDistributionEngine(ys, Z, list(xb = list(safeBackEnd)), type = "foo"), "'asymptotic' or 'exact'")
 
})
