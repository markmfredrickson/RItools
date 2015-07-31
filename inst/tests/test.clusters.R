################################################################################
# Tests for clustered treatment assignment
################################################################################

library("testthat")

context("Clustered treatment assignment")

test_that("Clusters must be aligned to treatment assignment and nested within strata", {
  set.seed(20130801)

  d <- data.frame(
      id          = 1:500,
      x           = rnorm(500),
      cluster     = rep(1:100, 5),
      strata.good = rep(1:5, 100),
      strata.bad  = sample(1:100, size = 500, replace = T),
      z.good      = rep(c(0,1,0,1,0), each = 100),
      z.bad       = sample(c(0,1), size = 500, replace = T))
     
  # check our assumptions that what is good and bad is really good and bad
  f <- function(x) { dim(x)[2] == 2 && all(apply(x, 1, function(i) { all(i != 0) })) }
  check.good <- by(d, INDICES = d$strata.good, FUN = function(i) { f(table(i$cluster, i$z.good)) })
  expect_true(all(check.good))

  check.bad.strata <- by(d, INDICES = d$strata.bad, FUN = function(i) { f(table(i$cluster, i$z.good)) })
  expect_true(any(!check.bad.strata))

  check.bad.z <- by(d, INDICES = d$strata.good, FUN = function(i) { f(table(i$cluster, i$z.bad)) })
  expect_true(any(!check.bad.z))

  # now make sure that xBalance complains appropriately

  # first, let's just check without clustering. Just look for non-errors
  # obviously, both good is fine
  xBalance(z.good ~ x + strata(strata.good), data = d)

  # more subtly, without the clustering, the "z.bad" is ok as it doesn't need to be aligned to clusters
  xBalance(z.bad ~ x + strata(strata.good), data = d)
  xBalance(z.good ~ x + strata(strata.bad), data = d)

  # now for our errors
  expect_error(xBalance(z.bad ~ x + strata(strata.good) + cluster(cluster), data = d), "cluster")
  expect_error(xBalance(z.good ~ x + strata(strata.bad) + cluster(cluster), data = d), "cluster")

  # finally, expect different results when clusters are included.
  resC <- xBalance(z.good ~ x + strata(strata.good) + cluster(cluster), data = d)
  resNC <- xBalance(z.good ~ x + strata(strata.good), data = d)

  expect_false(identical(resC, resNC))
  
  # when the cluster is a single observation, we should get exactly the same answers
  expect_equivalent(resNC, xBalance(z.good ~ x + strata(strata.good) + cluster(id), data = d))
})

