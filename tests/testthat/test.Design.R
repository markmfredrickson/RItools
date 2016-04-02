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

  d$'(weights)' = 1 # meet expectation of a weights column

  # checking input
  expect_error(makeDesign(~ x, data = d), "treatment")
  expect_error(makeDesign(z.bad ~ x + strata(strata.good) + cluster(cluster), data = d), "cluster")
  expect_error(makeDesign(z.good ~ x + strata(strata.bad) + cluster(cluster), data = d), "strata")
  expect_error(makeDesign(z.good ~ x + cluster(id) + cluster(cluster), data = d), "cluster")
  expect_error(makeDesign(z.good ~ x - 1, data = d, "stratification"))

  # actually testing that the output is as expected
  simple <- makeDesign(z.good ~ x, data = d)
  expect_equal(dim(simple@StrataFrame)[2], 1)
  expect_equivalent(simple@Covariates[, "x"], d$x)
  expect_equivalent(simple@Z, as.logical(d$z.good))
  expect_equal(nlevels(simple@Cluster), 500) # a cluster per individual
  
  clustered <- makeDesign(z.good ~ x + cluster(cluster), data = d)
  expect_equal(dim(clustered@StrataFrame)[2], 1)
  expect_true(nlevels(clustered@Cluster) > 1)

  clustStrata <- makeDesign(z.good ~ x + cluster(cluster) + strata(strata.good), data = d)
  expect_equal(dim(clustStrata@StrataFrame)[2], 2)

  # dropping the overall comparison
  expect_equal(dim(makeDesign(z.good ~ x + cluster(cluster) + strata(strata.good) - 1, data = d)@StrataFrame)[2], 1)
   
  ## More tests to write:
  # - All NA strata variables
  # - Missing z or cluster
  # - strata with extra levels but no observations (which can be safely dropped)
  #   (NB: extra levels tested upstream, in xBalance, as of commit 34861515; 
  #   see ./test.clusters.R  
})

test_that("Design to descriptive statistics", {
  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  d$'(weights)' = 1 # meet expectation of a weights column
  
  simple <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d)

  
  descriptives <- RItools:::designToDescriptives(simple)

  expect_is(descriptives, "array")
  expect_equal(dim(descriptives), c(4, 5, 2))

  # descriptives ignore clustering
  design.noclus <- RItools:::makeDesign(z ~ x + f + strata(s), data = d)
  expect_equal(descriptives, RItools:::designToDescriptives(design.noclus))
  
  # the strata should imply different stats
  expect_false(identical(descriptives[,,1], descriptives[,,2]))
  
  # ok, now checking that values are correct.
  expect_equal(mean(d$x[d$z == 1]), descriptives["x", "Treatment", "Unstrat"])
  expect_equal(mean(d$x[d$z == 0]), descriptives["x", "Control", "Unstrat"])

  # with equal sized strata, the the control/treatment means are the means of the the strata means
  expect_equal(mean(tapply(d$x[d$z == 1], d$s[d$z == 1], mean)), descriptives["x", "Treatment", "s"])
  expect_equal(mean(tapply(d$x[d$z == 0], d$s[d$z == 0], mean)), descriptives["x", "Control", "s"])

})

test_that("Issue 36: Descriptives with NAs, appropriate weighting", {
  ### Descriptives with missing covariates ###

  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      paired = rep(c(0,1), each = 250),
      z = rep(c(0,1), 250))
  d$'(weights)' = 1
  
  d.missing <- d

  d.missing$x[d.missing$x < -1] <- NA

  simple.all <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d)
  
  descriptives.all <- RItools:::designToDescriptives(simple.all)
  expect_equal(descriptives.all["x", "Treatment", "Unstrat"], mean(d$x[d$z == 1]))
  expect_equal(descriptives.all["x", "Treatment", "s"], mean(d$x[d$z == 1 & !is.na(d$s)]))

  simple.missing <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d.missing)
  
  descriptives.missing <- RItools:::designToDescriptives(simple.missing)

  with(d.missing,
       expect_equal(descriptives.missing["x", "Treatment", "Unstrat"],
                    mean(x[z == 1], na.rm = TRUE)))

  with(d.missing,
       expect_equal(descriptives.missing["x", "Treatment", "s"],
                    mean(x[z == 1 & !is.na(s)], na.rm = TRUE)))

  expect_false(identical(descriptives.all, descriptives.missing))

  # ETT weighting
  design.paired   <- RItools:::makeDesign(z ~ x + f + strata(paired) + strata(s), data = d) 
  descriptives.paired <- RItools:::designToDescriptives(design.paired)

  with(d, expect_equal(descriptives.paired["x", "Control", "paired"], mean(x[z == 0])))
  with(d, expect_equal(descriptives.paired["x", "Control", "Unstrat"], mean(x[z == 0])))
  with(d, expect_false(identical(descriptives.paired["x", "Control", "s"], mean(x[z == 0]))))
})

test_that("Aggegating designs by clusters", {
  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  d$'(weights)' <- 1
  # grab a bunch of rows and duplicate them
  d <- rbind(d, d[sample(1:dim(d)[1], size = 100), ])
  
  # swapping around the data to make sure order doesn't mask bugs
  d <- d[sample(1:dim(d)[1]), ]


  design <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d)
  aggDesign <- RItools:::aggregateDesign(design) 

  # one row per cluster, with columns x, fa, fb, fc
  expect_equal(dim(aggDesign@Covariates), c(100, 4))

  # now spot check some cluster totals of totals
  expect_equal(aggDesign@Covariates[1, ], colMeans(design@Covariates[design@Cluster == 1,]))

})

test_that("aggregateDesign treats NA covariates as 0's" ,{
  dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x=rep(c(NA, 1), c(3,2))
                    )
  dat$'(weights)' <- 1

  design <- RItools:::makeDesign(z~x+strata(strat)+cluster(clus)-1, dat)
  aggDesign <- RItools:::aggregateDesign(design)
  expect_equal(aggDesign@Covariates[,'x'], c(0, 0, 1, 1) )
  expect_equal(aggDesign@Eweights[,colnames(aggDesign@Covariates)=='x'], c(0,0,1,1) )
})

test_that("NA flags are optional", {
  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      paired = rep(c(0,1), each = 250),
      z = rep(c(0,1), 250))
  d$'(weights)' <- 1

  d$x[sample.int(500, size = 10)] <- NA

  design.flags   <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d)
  design.noFlags <- RItools:::makeDesign(z ~ x + f + strata(s), data = d, include.NA.flags = FALSE)

  expect_equal(dim(design.flags@Covariates)[2], 5)
  expect_equal(dim(design.noFlags@Covariates)[2], 4)
})

test_that("Aggregation of element weights to cluster level",{
  ##set.seed(20130801)

  d.short <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  
  ## grab a bunch of rows and duplicate them
  newrows <- c(1:nrow(d.short), sample(1:nrow(d.short), size=100, replace=T) )
  d.tall <- d.short[newrows,]
  d.tall$'(weights)' <- 1
  d.short$'(weights)' <- as.vector(table(newrows))

  design.tall <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d.tall)
  aggDesign.tall <- RItools:::aggregateDesign(design.tall) 

    design.short <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d.short)
  aggDesign.short <- RItools:::aggregateDesign(design.short) 

  ## we should wind up in the same place.
  expect_equal(aggDesign.tall, aggDesign.short)

  d2 <- d.tall
  d2$'(weights)' <- 2
  design.d2 <- RItools:::makeDesign(z ~ x + f + strata(s) + cluster(c), data = d2)
  aggDesign.d2 <- RItools:::aggregateDesign(design.d2)
  
  expect_equal(2* aggDesign.tall@Eweights, aggDesign.d2@Eweights)
  expect_equal(aggDesign.tall@Covariates, aggDesign.d2@Covariates) 
          })

