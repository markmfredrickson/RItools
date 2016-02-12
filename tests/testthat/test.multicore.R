################################################################################
# Parallel package and mclapply tests
################################################################################

library(testthat)

context("Mulicore package use")

test_that("Check parallel on/off", {

  if ("parallel" %in% loadedNamespaces()) {
    unloadNamespace("parallel")  
  }

  expect_false(parallelLoaded())

  serial <- getApplyFunction(1:10^2, sqrt)

  library(parallel)

  expect_true(parallelLoaded())

  parallel <- getApplyFunction(1:10^2, sqrt)

  expect_equal(serial, parallel)
  
})

test_that("Runs work in all systems", {
  
  if ("parallel" %in% loadedNamespaces()) {
    unloadNamespace("parallel")  
  }

  expect_false(parallelLoaded())
  
  set.seed(20110620)
  tau <- 10
  n <- 8
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  set.seed(123456)
  time.nm <- system.time(res.nomulti <- RItest(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(8, 12, .25))))

  library(parallel)
  
  set.seed(123456)
  time.mul <- system.time(res.multi <- RItest(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(8, 12, .25))))

  expect_equal(res.nomulti, res.multi)
  
  # This fails about 1/2 the time on a 2 core box
  # it is not a requirement that mcapply be faster (though it would be nice)
  # expect_true(time.mul[[3]] < time.nm[[3]]) # elapsed time

  opts <- options("RItools-apply" = mclapply)
  
  set.seed(123456)
  res.explicit <- RItest(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(8, 12, .25)))
  
  expect_equal(res.nomulti, res.explicit)

  options(opts)
})

# I am unable to change "locked bindings" of mclapply or lapply
# so commenting this out for now.
### test_that("mcapply called", {
###   if ("parallel" %in% loadedNamespaces()) {
###     unloadNamespace("parallel")  
###   }
### 
###   expect_false(parallelLoaded())
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
###   library(parallel)
###   old.mclapply <- mclapply
###   called <- FALSE
###   mclapply <<- function(...) { called <<- T ; error("Good!") }
### 
### 
###   expect_error(RItest(R, Z, mann.whitney.u))
###   expect_true(called)
###   
###   mclapply <<- old.mclapply
### })
