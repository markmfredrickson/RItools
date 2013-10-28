################################################################################
# Tests for functions in src/algorithms.cpp
################################################################################

library("testthat")

context("C++ algorithms")

test_that("running models of different param numbers", {
  set.seed(20131026)
  n <- 1000
  y <- rnorm(n) 

  z <- rep(c(1,0), n/2)

  zs <- as.matrix(simpleRandomSampler(n, sum(z))(100)$samples)
  
  p1 <- -5:5
  res1 <- runModel1(meanDifference, y, z, 
                   zs, constant.additive.model, p1) 

  expect_equal(length(res1), length(p1))

  twoparamMod <- function(y, z, a, b) {
    y - z * (a + b)
  }
  p2 <- c(100, 200, 300)

  res2 <- runModel2(meanDifference, y, z,
                    zs, twoparamMod, p1, p2)

  expect_equal(length(res2), length(p1) * length(p2))


  ### with the basics out of the way, lets check correctness by implementing the same technique purely in R
  
  test.stat <- meanDifference
  # getting the p-value
  doit <- function(m) {
    function(...) { 
      # first, for each model, adjust the observed y with the observed
      # z
      adjusted.y <- m(y, z, ...)
      obs.t <- test.stat(adjusted.y, z)

      # check if there is a backend function for this test.statistic
      # if (type == "asymptotic") {
      #     if (inherits(test.stat, "AsymptoticTestStatistic")) {
      #       # basically, pass the adjusted y and everything else
      #       # to the backend and let it return a RandomizationDistribution
      #       # or subclass (probably a good idea)
      #       # backends may not honor samples, p.value, or summaries
      #       return(test.stat@asymptotic(adjusted.y, z))
      #     }

      #   stop("No asymptotic backend exists for the test statistic")
      # }

      # no backend, so use the standard approach
      # Possible todo: pull this out into its own "backend" function

      # now iterate over the randomizations, using the adjusted y
      # if we can, we'll try to call the test stat at the C++ level using a function pointer
      tsf <- test.stat
      if (inherits(test.stat, "RcppFunction")) {
        tsf <- test.stat@RcppFn
      }

      this.p <- computeTestStatPval(tsf, adjusted.y, obs.t, zs)
        
      return(this.p)
    }
  }

  res1.farray <- farray(doit(constant.additive.model), list(p1))
  res2.farray <- farray(doit(twoparamMod), list(p1, p2))

  expect_equal(res1, as.vector(res1.farray))
  expect_equal(res2, as.vector(res2.farray))
})
