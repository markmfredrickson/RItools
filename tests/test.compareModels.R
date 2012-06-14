################################################################################
# Tests of the compareModels function
################################################################################  

library(testthat)

context("compareModels function")

test_that("testing a single model", {
  # This is a smaller version of the lady tasting tea
  # 20 * 20 is 400 total loops, versus 70 * 70 = 4900 loops

  number.correct <- function(guesses, z) { sum(z == guesses) / 2 }
  
  sampler <- simpleRandomSampler(treated = 3, total = 6)

  res <- compareModels(models = list(sharp = sharp.null.model),
                       repetitions = 20,
                       samples = 20,
                       test.statistic = number.correct,
                       sampler = sampler,
                       uniformity = c(1,0,1,0,0,1),
                       p.value = upper.p.value)

  expect_equal(class(res), "array") 
  expect_equal(dim(res), c(1,1,20))

  # collect the p-values
  pvs <- res[1,1,]
  expect_equal(as.numeric(table(pvs)), c(1, 9, 9, 1))

  expect_equal(names(pvs), as.character(1:20))
  
})
