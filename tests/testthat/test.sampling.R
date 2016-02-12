################################################################################
# Sampling mechanisms
################################################################################

library(testthat)

context("Sampling")

test_that("simpleRandomSample for fixed number of treated within blocks", {
  set.seed(20131010)
  sampler <- simpleRandomSampler(total = 8, treated = 4) # 70 possible
  
  expect_equal(dim(sampler(100)$samples)[2], 70) # do not return more than 70
  expect_equal(dim(sampler(4)$samples)[2], 4) # less than total
  expect_false(all(sampler(4)$samples == sampler(4)$samples)) # randomize

  # next up we have (4 choose 2) * (4 choose 1) = 24 possible z's
  blocked <- simpleRandomSampler(total = c(4,4), treated = c(2,1))
  expect_equal(dim(blocked(100)$samples)[2], 24)
  expect_equal(dim(blocked(10)$samples)[2], 10)

  # data should be in block order 11112222
  tmp <- blocked(25) # enumerate
  expect_true(all(colSums(tmp$samples[1:4,]) == 2))
  expect_true(all(colSums(tmp$samples[5:8,]) == 1))

  # alternative invocation takes z and b arguments
  # can be useful when data are not in block order
  b <- c(1,2,3,1,2,3,3,2,1,3,2,1) # 3 blocks, 4 units each
  z <- c(1,1,1,0,0,0,0,0,1,0,0,0) # 2 in block 1, 1 in block 2 and 3
  indexed <- simpleRandomSampler(z = z, b = b)
  indexed.res <- indexed(100)$samples
  expect_equal(dim(indexed.res)[2], 6 * 4 * 4)
  
  # check that the randomizations are in data order, not block order
  expect_true(all(colSums(indexed.res[c(1,4,9,12),]) == 2))
  expect_true(all(colSums(indexed.res[c(2,5,8,11),]) == 1))

  # error handling
  expect_error(simpleRandomSampler(total = c(2,3), treated = c(1,3)))
  expect_error(simpleRandomSampler(total = c(2,3), treated = c(1,0)))
  expect_error(simpleRandomSampler(total = c(2,3,4), treated = c(1,2)))

  # one of total/b and treated/z must be passed (they can be cross inferred)
  # valid
  simpleRandomSampler(total = c(2,3), z = c(1,0,0,0,1))
  simpleRandomSampler(b = c(1,1,2,2,1), treated = c(2, 1))
  # invalid -- don't look terrible, but let's not allow them
  # the first doesn't have enough info, the second could be perhaps inferred
  # but let's not allow it either
  expect_error(simpleRandomSampler(total = c(2,3), b = c(1,1,2,2,2)))
  expect_error(simpleRandomSampler(z = c(1,0,0,1,0,0), treated = c(1,1)))

  # blocks as factors
  bf <- simpleRandomSampler(z = c(1,0,0,1), b = as.factor(c("a", "b", "a", "b")))
  expect_equal(dim(bf(100)$samples), c(4, 4))
  expect_true(all(colSums(bf(100)$samples[c(1,3),]) == 1))

  # NAs are an error. The user is responsible for cleaning NAs before sampling.
  expect_error(simpleRandomSampler(z = c(1,0,0,1,0), b = c(1,2,1,2,NA)), "NAs")
})

test_that("multinomialSampler for independent bernoulli draws", {

  set.seed(20131010)
  sampler <- independentProbabilitySampler(4) 

  # max 2^4 = 16 randomizations (including all zero and all one)
  expect_equal(dim(sampler(100)$samples)[2], 16)
  expect_equal(dim(sampler(4)$samples)[2], 4)
  expect_equal(sum(sampler(100)$samples), 32)
  expect_true(all(sampler(100)$weight == 1))
  
  unequal.sampler <- independentProbabilitySampler(2, c(0.25, 0.5))
  expect_true(all(unequal.sampler(10)$weight %in% c(0.75 * 0.5, 0.25 * 0.5))) 

})
