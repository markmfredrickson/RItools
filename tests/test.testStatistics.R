################################################################################
# Testing the behind the scenes engine 
################################################################################

library(testthat)

context("Test statistics, functions and objects")

# Turning off until blocking gets back in the game
#test_that("Harmonic Mean Difference", {
#  
#  set.seed(20111018)
#  n <- 200 # pick a bigger n than any of our other tests, as this should be fast
#  Z <- rep(c(0,1), n/2)
#  B <- rep(1:4, each =  n/4)
#  ys <- rnorm(n) + Z + B/8 # small block effect, larger treatment effect
#  
#  hmd <- xBalance(Z ~ ys, data = data.frame(Z, ys), strata = as.factor(B), report =
#    "adj.mean.diffs")$results[,1,1]
#
#  expect_equal(harmonic.mean.difference(ys, Z, B), hmd)
#
#  # the model, a location shift (does not fit reality as it ignores block
#  # effects, but that's fine. no model is perfect.)
#  tau1 <- function(y, z, b) { modelOfEffect(constant.additive.model, ys, z, b, tau = 1)}
#  
#  # samples should be ignored, set low to keep the test short if there is an error
#  res.xb <- randomizationDistributionEngine(ys, Z, 
#                                            list(xb = list(harmonic.mean.difference, tau1)), 
#                                            blocks = B,
#                                            samples = 1,
#                                            summaries = "z.scores",
#                                            type = "asymptotic") 
#  dst <- res.xb$xb
#
#  expect_equal(dim(dst), c(2,3))
#  expect_equal(colnames(dst), c("statistic", "p.value", "z"))
#
#  # now use the exact version and check our results
#  # by default type = "exact"
#  res.exact <- randomizationDistributionEngine(ys, Z, 
#                                            list(xb = list(harmonic.mean.difference, tau1)), 
#                                            blocks = B,
#                                            samples = 1)
#
#  dst.x <- res.exact$xb
#
#  expect_equal(dst$statistic, dst.x$statistic)
#
#})

test_that("Wilcox.test", {
  
  set.seed(20111018)
  n <- 200 # pick a bigger n than any of our other tests, as this should be fast
  Z <- rep(c(0,1), n/2)
  B <- rep(1:4, each =  n/4)
  ys <- rnorm(n) + Z + B/8 # small block effect, larger treatment effecto
  tau1 <- function(y, z) { constant.additive.model(ys, z, tau = 1) }
  
  res.wt <- randomizationDistributionEngine(ys, 
                                            Z, 
                                            list(wt = list(mann.whitney.u, tau1)),  
                                            samples = 1,
                                            type = "asymptotic")

  dst <- res.wt$wt

  expect_equal(dim(dst), c(2,2))
  expect_equal(colnames(dst), c("statistic", "p.value"))

  # the exact version
  res.exact <- randomizationDistributionEngine(ys, 
                                            Z, 
                                            list(wt = list(mann.whitney.u, tau1)),  
                                            samples = 1)

  dst.x <- res.exact$wt
  expect_equal(dst$statistic, dst.x$statistic)

})

test_that("Quintile Differences", {

  y <- c(1:5, 2:6)
  z <- c(rep(0,5), rep(1,5))

  median.difference <- quantileAbsoluteDifference(0.5)
  expect_equal(median.difference(y, z), 1)

  y <- c(1,2,3,4,5, 1,2,3.5,4,5)
  quartile.difference <- quantileAbsoluteDifference(c(0, 0.25, 0.5, 0.75, 1))
  expect_equal(quartile.difference(y, z), 0.5)

})

test_that("Subgroup Analysis", {
  y <- c(1:5, 2:6)
  z <- c(rep(0,5), rep(1,5))
  x <- c(rep(1,8), 0,0)

  xdiff <- subsetStatistic(mean.difference, x == 1)
  
  expect_equal(xdiff(y, z), 0)
})

test_that("KS Test Statistic", {
  
  y <- c(seq(1,11, by = 2), seq(10,20, by = 2))
  z <- c(rep(0,6), rep(1,6))
  
  ### Agree with an exact computation from ks.test
  res.ks <- ks.test(y[z == 1], y[z == 0])

  expect_equal(as.numeric(res.ks$statistic), ksTestStatistic(y, z))

  res.ksfn <- RItest(y, z, ksTestStatistic, p.value = upper.p.value)

  expect_equivalent(res.ks$statistic, res.ksfn$statistic)
  expect_equal(res.ks$p.value, res.ksfn$p.value)

  ### Next tests making sure that we can flip things around and get the right values
  res.ks.flip <- ks.test(y[z == 0], y[z == 1])

  expect_equal(as.numeric(res.ks.flip$statistic), ksTestStatistic(y, rev(z)))

  res.ksfn.flip <- RItest(y, rev(z), ksTestStatistic, p.value = upper.p.value)
  

  expect_equivalent(res.ks.flip$statistic, res.ksfn.flip$statistic)
  expect_equal(res.ks.flip$p.value, res.ksfn.flip$p.value)

  ### Compare the backend results to calling ks.test explicity for large data. Should be the same

  set.seed(20120423)
  y0 <- rnorm(10000)
  z <- rep(c(0,1), 5000)
  y <- y0 + 0.05 * z

  res.ks.asym <- ks.test(y[z == 1], y[z == 0], exact = F) # be sure to use asymptotics
  res.ksfn.asym <- RItest(y, z, ksTestStatistic, type = "asymptotic")
  
  expect_equivalent(res.ks.asym$statistic, res.ksfn.asym$statistic)
  expect_equal(res.ks.asym$p.value, res.ksfn.asym$p.value)


})
