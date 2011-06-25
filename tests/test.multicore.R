################################################################################
# Multicore/mclapply tests
################################################################################

library(testthat)

context("Mulicore package use")

test_that("Check multicore on/off", {

  if ("multicore" %in% loadedNamespaces()) {
    unloadNamespace("multicore")  
  }

  expect_false(multicoreLoaded())

  expect_true(is.null(options("RItools-apply")[[1]]))
  expect_equal(getApplyFunction(), lapply)

  library(multicore)

  expect_true(multicoreLoaded())
  expect_equal(getApplyFunction(), mclapply)

  opts <- options("RItools-apply" = sum) # really bad choice :-)
  expect_equal(getApplyFunction(), sum)
  options(opts)
  
})

test_that("Runs work in all systems", {
  
  if ("multicore" %in% loadedNamespaces()) {
    unloadNamespace("multicore")  
  }

  expect_false(multicoreLoaded())
  
  set.seed(20110620)
  tau <- 10
  n <- 8
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  set.seed(123456)
  time.nm <- system.time(res.nomulti <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(8, 12, .25))))

  library(multicore)
  
  set.seed(123456)
  time.mul <- system.time(res.multi <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(8, 12, .25))))

  expect_equal(res.nomulti, res.multi)
  expect_true(time.mul[[3]] < time.nm[[3]]) # elapsed time

  opts <- options("RItools-apply" = mclapply)
  
  set.seed(123456)
  res.explicit <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(8, 12, .25)))
  
  expect_equal(res.nomulti, res.explicit)

  options(opts)
})

# I am unable to change "locked bindings" of mclapply or lapply
# so commenting this out for now.
### test_that("mcapply called", {
###   if ("multicore" %in% loadedNamespaces()) {
###     unloadNamespace("multicore")  
###   }
### 
###   expect_false(multicoreLoaded())
###   expect_true(is.null(options("RItools-apply")[[1]]))
### 
###   set.seed(20110620)
###   tau <- 10
###   n <- 6
###   Yc <- rnorm(n)
###   Yt <- Yc + tau
###   Z <- rep(0, n)
###   Z[sample.int(n, n/2)] <- 1
###   R <- Z * Yt + (1 - Z) * Yc
### 
###   library(multicore)
###   old.mclapply <- mclapply
###   called <- FALSE
###   mclapply <<- function(...) { called <<- T ; error("Good!") }
### 
### 
###   expect_error(parameterizedRandomizationDistribution(R, Z, mann.whitney.u))
###   expect_true(called)
###   
###   mclapply <<- old.mclapply
### })
