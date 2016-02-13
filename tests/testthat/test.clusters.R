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
      z.good      = rep(c(0,1), 250),
      z.bad       = sample(c(0,1), size = 500, replace = T))
     
  # check our assumptions that what is good and bad is really good and bad
  check <- function(tbl) {
    cluster <- apply(tbl, 2, function(i) { sum(i != 0) == 1 })
    strata  <- apply(tbl, 3, function(i) { dim(i)[1] == 2 && all(rowSums(i) > 0) })
    return(all(cluster) && all(strata))
  }

  expect_true(check(table(d[, c("z.good", "cluster", "strata.good")])))
  expect_false(check(table(d[, c("z.bad", "cluster", "strata.good")])))
  expect_false(check(table(d[, c("z.good", "cluster", "strata.bad")])))
  expect_true(check(table(d[, c("z.good", "id", "strata.good")])))
  

  # now make sure that xBalance complains appropriately

  # first, let's just check without clustering. Just look for non-errors
  # obviously, both good is fine
  xBalance(z.good ~ x + strata(strata.good), data = d)

  # more subtly, without the clustering, the "z.bad" is ok as it doesn't need to be aligned to clusters
  xBalance(z.bad ~ x + strata(strata.good), data = d)

  # now for our errors
  expect_error(xBalance(z.good ~ x + strata(strata.bad), data = d), "strata")
  expect_error(xBalance(z.bad ~ x + strata(strata.good) + cluster(cluster), data = d), "cluster")
  expect_error(xBalance(z.good ~ x + strata(strata.bad) + cluster(cluster), data = d), "cluster")

  # expect different results when clusters are included.
  resC <- xBalance(z.good ~ x + strata(strata.good) + cluster(cluster), data = d)
  resNC <- xBalance(z.good ~ x + strata(strata.good), data = d)

  expect_false(identical(resC, resNC))
  
  # when the cluster is a single observation, we should get exactly the same answers
  expect_equivalent(resNC, xBalance(z.good ~ x + strata(strata.good) + cluster(id), data = d))

  # only one cluster arg allowed
  expect_error(xBalance(z.good ~ x + cluster(cluster) + cluster(id), data = id))
})

