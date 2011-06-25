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
