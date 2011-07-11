################################################################################
# Testing the parameterizedRandomizationDistribution against wilcox.test
################################################################################

library(testthat)

context("Wilcox-like results")

##Unpaired
test_that("Get same point estimate", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  expect_equal(sum(Z), 6)

  res.wilcox <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T)

  res.prd <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = c(-10, 9, 10)))

  # the wilcox stat is named "W", but is otherwise the same
  # the result of the test stat applied to the observed data it the first
  # entry in the sharp null object
  expect_equal(res.prd[[1,1]], res.wilcox$statistic[[1]])

  s <- summary(res.prd)
  expect_equal(res.wilcox$p.value, s@sharp.null.p)
  
  res.wilcox.m10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = -10)
  res.wilcox.10 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 10)
  res.wilcox.9 <- wilcox.test(R[Z == 1], R[Z == 0], exact = T, conf.int = T, mu = 9)

  pvs <- na.omit(p.values(res.prd)) # na.omit drops the sharp null, which complicates
  
  expect_equal(pvs[pvs$tau == 10, "p"], res.wilcox.10$p.value)
  expect_equal(pvs[pvs$tau == -10, "p"], res.wilcox.m10$p.value)
  expect_equal(pvs[pvs$tau == 9, "p"], res.wilcox.9$p.value)

  # zeroing in on the point estimate
  step <- 0.01
  res.prd.intval <- parameterizedRandomizationDistribution(R, Z, mann.whitney.u, 
    constant.additive.model, list(tau = seq(9,10.5, step)))

  # this problem does not have a unique point estimate. Typical HL would then take the
  # mid point of the interval. In our case, we'll take the mean and call it close enough.
  s.intval <- summary(res.prd.intval)
  thept <- mean(s.intval@point.estimate$tau)
  expect_true(abs(thept - res.wilcox$estimate[[1]]) < step) # [[]] to remove name
})

test_that("Get same point estimate for Paired Data (where the pairing is not predictive of outcomes but more or less random and mainly exists to constrain the possible randomizations)", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  nB <- n/2 ## number of blocks, here, pairs
  
  Yc <- rnorm(n)
  Yt <- Yc + tau
  Z <- rep(c(0,1), n/2)
  ##Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc
  B <- gl(nB,2)
  
  expect_equal(sum(Z), 6)


align.by.block<-function (x, block, fn = mean) {
    ##"aligns" values of a variable within blocks/strata/sets.
    ##Default is to subtract off the mean.
    unsplit(lapply(split(x, block), function(x) {x - fn(x)}), block)
  }


paired.constant.additive.model<-function(ys,z,b,tau){ ##where z is really V_b=Z_1b - Z_2b
  (ys*z)-(tau*z)
}
  
  ##Paired version
  R.paired <-sapply(split(data.frame(R=R,Z=Z),B),function(dat){
    with(dat,R[Z==1]-R[Z==0])
  })
  Z.paired<-sapply(split(Z,B),function(zb){zb[2]-zb[1]})

  ##Notice that this is what happens inside of wilcox.test.default around line 23
  ##if (paired) {
  ##          if (length(x) != length(y)) 
  ##              stop("'x' and 'y' must have the same length")
  ##          OK <- complete.cases(x, y)
  ##          x <- x[OK] - y[OK]
  ##          y <- NULL
  ##      }

  expect_equal(mann.whitney.u(R.paired,rep(1,n/2)),
               paired.sgnrank.sum(R,Z,B))

  ##For comparison (using means)
  lm1<-lm(R~Z+B)
  lm2<-lm(R.aligned~Z)
  lm3<-lm(R.paired~1)
  expect_equal(coef(lm1)["Z"],coef(lm2)["Z"],coef(lm3),check.attributes=FALSE)

  ##and same stat inf
  confint(update(lm2,.~.+B),parm="Z")
  confint(lm1,parm="Z")
  confint(lm3)

  
  ##Aligned Version (we'd have to change the model of effects to reflect the paired nature of the data)
  R.aligned <- align.by.block(R,B)
  
  res.wilcox <- wilcox.test(R.paired, exact = T, conf.int = T)
  ##notice that because of sorting by B and Z res.wilcox is same as wilcox.test(R[Z==1],R[Z==0],exact=T,conf.int=T,paired=TRUE),

  res.prd <- parameterizedRandomizationDistribution(R.aligned, Z, paired.sgnrank.sum, #mann.whitney.u, 
    constant.additive.model, list(tau = c(-10, 9, 10))) ##works for both R and R.aligned

  # the wilcox stat is named "W", but is otherwise the same
  # the result of the test stat applied to the observed data it the first
  # entry in the sharp null object
  expect_equal(res.prd[[1,1]], res.wilcox$statistic[[1]])

  s <- summary(res.prd)
  expect_equal(res.wilcox$p.value, s@sharp.null.p)
  
  res.wilcox.m10 <- wilcox.test(R.paired, exact = T, conf.int = T, mu = -10)
  res.wilcox.10 <-  wilcox.test(R.paired, exact = T, conf.int = T, mu = 10)
  res.wilcox.9 <-   wilcox.test(R.paired, exact = T, conf.int = T, mu = 9)

  pvs <- na.omit(p.values(res.prd)) # na.omit drops the sharp null, which complicates
  
  expect_equal(pvs[pvs$tau == 10, "p"], res.wilcox.10$p.value)
  expect_equal(pvs[pvs$tau == -10, "p"], res.wilcox.m10$p.value)
  expect_equal(pvs[pvs$tau == 9, "p"], res.wilcox.9$p.value)

  # zeroing in on the point estimate
  step <- 0.01
  res.prd.intval <- parameterizedRandomizationDistribution(R.aligned, Z, paired.sgnrank.sum, ##mann.whitney.u, 
    constant.additive.model, list(tau = seq(9,10.5, step)))

  # this problem does not have a unique point estimate. Typical HL would then take the
  # mid point of the interval. In our case, we'll take the mean and call it close enough.
  s.intval <- summary(res.prd.intval)
  thept <- mean(s.intval@point.estimate$tau)
  expect_true(abs(thept - res.wilcox$estimate[[1]]) < step) # [[]] to remove name
})
