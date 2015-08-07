################################################################################
# Tests for Design objects
################################################################################

library("testthat")

context("Design objects")

test_that("Creating design objects", {
  set.seed(20130801)

  d <- data.frame(
      id          = 1:500,
      x           = rnorm(500),
      cluster     = rep(1:100, 5),
      strata.good = rep(c(1:4, NA), 100),
      strata.bad  = sample(1:100, size = 500, replace = T),
      z.good      = rep(c(0,1), 250),
      z.bad       = sample(c(0,1), size = 500, replace = T))

  # checking input
  expect_error(makeDesign(~ x, data = d), "treatment")
  expect_error(makeDesign(z.bad ~ x + strata(strata.good) + cluster(cluster), data = d), "cluster")
  expect_error(makeDesign(z.good ~ x + strata(strata.bad) + cluster(cluster), data = d), "strata")
  expect_error(makeDesign(z.good ~ x + cluster(id) + cluster(cluster), data = d), "cluster")
  expect_error(makeDesign(z.good ~ x - 1, data = d, "stratification"))

  # actually testing that the output is as expected
  simple <- makeDesign(z.good ~ x, data = d)
  expect_equal(dim(simple@Strata)[2], 1)
  expect_equivalent(simple@Covariates[, "x"], d$x)
  expect_equivalent(simple@Z, as.logical(d$z.good))
  expect_equal(nlevels(simple@Cluster), 1)
  
  clustered <- makeDesign(z.good ~ x + cluster(cluster), data = d)
  expect_equal(dim(clustered@Strata)[2], 1)
  expect_true(nlevels(clustered@Cluster) > 1)

  clustStrata <- makeDesign(z.good ~ x + cluster(cluster) + strata(strata.good), data = d)
  expect_equal(dim(clustStrata@Strata)[2], 2)

  # dropping the overall comparison
  expect_equal(dim(makeDesign(z.good ~ x + cluster(cluster) + strata(strata.good) - 1, data = d)@Strata)[2], 1)
   
  ## More tests to write:
  # - All NA strata variables
  # - Missing z or cluster
  # - strata with extra levels but no observations (which can be safely dropped)
})

test_that("Design to descriptive statistics", {
  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  design <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d)
  simple <- RItools:::weightedDesign(design)
  
  descriptives <- RItools:::weightedDesignToDescriptives(simple)

  expect_is(descriptives, "array")
  expect_equal(dim(descriptives), c(4, 5, 2))

  # descriptives ignore clustering
  design.noclus <- RItools:::makeDesign(z ~ x + f + strata(s), data = d)
  expect_equal(descriptives, weightedDesignToDescriptives(weightedDesign(design.noclus)))
  
  # the strata should imply different stats
  expect_false(identical(descriptives[,,1], descriptives[,,2]))
  
  # ok, now checking that values are correct.
  expect_equal(mean(d$x[d$z == 1]), descriptives["x", "Treatment", "Unstrat"])
  expect_equal(mean(d$x[d$z == 0]), descriptives["x", "Control", "Unstrat"])

  # with equal sized strata, the the control/treatment means are the means of the the strata means
  expect_equal(mean(tapply(d$x[d$z == 1], d$s[d$z == 1], mean)), descriptives["x", "Treatment", "s"])
  expect_equal(mean(tapply(d$x[d$z == 0], d$s[d$z == 0], mean)), descriptives["x", "Control", "s"])
})

test_that("Aggegating designs by clusters", {
  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  # grab a bunch of rows and duplicate them
  d <- rbind(d, d[sample(1:dim(d)[1], size = 100), ])
  
  # swapping around the data to make sure order doesn't mask bugs
  d <- d[sample(1:dim(d)[1]), ]


  design <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d)
  aggDesign <- RItools:::aggregateDesign(design) 

  # one row per cluster, with columns x, fa, fb, fc, and cluster.size
  expect_equal(dim(aggDesign@Covariates), c(100, 5))

  # properly aggregating up cluster counts
  expect_equivalent(aggDesign@N, as.vector(table(d$c))[aggDesign@Cluster])

  # now spot check some cluster totals of totals
  expect_equal(aggDesign@Covariates[1, -1], colSums(design@Covariates[design@Cluster == aggDesign@Cluster[1],]))

})
