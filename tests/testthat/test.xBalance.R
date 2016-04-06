################################################################################
# Tests for xBalance function
################################################################################
## if working interactively in inst/tests you'll need
## library(RItools, lib.loc = '../../.local')
library("testthat")

context("xBalance Functions")

test_that("xBal univariate desriptive means agree w/ lm",{
    set.seed(20160406)
    n <- 7 
     dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
     dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )

     lm1 <- lm(x1~z, data=dat)
     xb1 <- xBalance(z~x1+strata(s), data=dat, report=c("adj.mean.diffs"))
     expect_equal(xb1$results["x1", "adj.diff", "Unstrat"], coef(lm1)["z"], check.attributes=F)

     ## try to match default ETT weighting
     pihat <- fitted(lm(z~s, data=dat))     
     lm2a <- lm(x1~z+s, data=dat, weights=ifelse(pihat==1,1, (1-pihat)^-1))

     expect_equal(xb1$results["x1", "adj.diff", "s"], coef(lm2a)["z"], check.attributes=F)


})

test_that("xBal univariate inferentials agree w/ conditional logistic Rao score test",{
    library(survival)
    set.seed(20160406)
    n <- 7 
     dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
     dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )
    xb1b <- xBalance(z~x1+strata(s), data=dat, report=c("z.scores"))
     cl1 <- clogit(z~x1, data=dat)
     cl2 <- clogit(z~x1+strata(s), data=dat)


    expect_equal(summary(cl1)$sctest['test'],(xb1b$results["x1", "z", "Unstrat"])^2 , check.attributes=F)
    expect_equal(summary(cl2)$sctest['test'],(xb1b$results["x1", "z", "s"])^2 , check.attributes=F)
}
          )

test_that("xBalance returns covariance of tests", {
  set.seed(20130801)
  n <- 500

  library(MASS)
  xs <- mvrnorm(n,
                mu = c(1,2,3),
                Sigma = matrix(c(1, 0.5, 0.2,
                    0.5, 1, 0,
                    0.2, 0, 1), nrow = 3, byrow = T))

  p <- plogis(xs[,1]- 0.25 * xs[,2] - 1)
  z <- rbinom(n, p = p, size = 1)
  s <- rep(c(0,1), each = n/2)

  dat <- cbind(z, xs, s)


  # we use ETT weighting here to correspond to the weighting scheme used
  # in the descriptives section
  res <- xBalance(z ~ . + strata(s),
                  data = as.data.frame(dat),
                  stratum.weights = RItools:::effectOfTreatmentOnTreated,
                  report = 'all')

  tcov <- attr(res$overall, "tcov")

  expect_false(is.null(tcov))

  expect_equal(length(tcov), 2)
  expect_equal(dim(tcov[[1]]), c(4,4))

})

test_that("Passing post.alignment.transform, #26", {
  data(nuclearplants)

  # Identity shouldn't have an effect
  res1 <- xBalance(pr ~ ., data=nuclearplants)
  res2 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = function(x) x)

  expect_true(all.equal(res1, res2)) ## allow for small numerical differences

  res3 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = rank)

  expect_true(all(dim(res1$results) == dim(res3$results)))

  expect_error(xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = mean),
               "Invalid post.alignment.transform given")

  res4 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = rank, report="all")
  res5 <- xBalance(pr ~ ., data=nuclearplants, report="all")

  expect_false(isTRUE(all.equal(res4,res5)))

  # a wilcoxon rank sum test, asymptotic and w/o continuity correction
  res6 <- xBalance(pr ~ cost, data=nuclearplants, post.alignment.transform = rank, report="all")

  expect_equal(res6$results["cost", "p", "Unstrat"],
               wilcox.test(cost~pr, data=nuclearplants, exact=FALSE, correct=FALSE)$p.value)

  # w/ one variable, chisquare p value should be same as p value on that variable
  expect_equal(res6$results["cost", "p", "Unstrat"],
               res6$overall["Unstrat","p.value"])

  # to dos: test combo of a transform with non-default stratum weights.

})

test_that("NA in stratify factor are dropped", {
  data(nuclearplants)

  n2 <- nuclearplants
  n2 <- rbind(n2, n2[1,])
  n2$pt[1] <- NA

  f <- function(d) {
    xBalance(pr ~ . - pt + strata(pt) - 1, data = d)
  }

  xb1 <- f(nuclearplants)
  xb2 <- f(n2)

  expect_equal(xb1, xb2)
})

test_that("Use of subset argument", {
  data(nuclearplants)

  n2 <- nuclearplants
  n2 <- rbind(n2, n2[1,])
  n2[nrow(nuclearplants)+1, "pt"] <- 2

  xb1 <- xBalance(pr ~ . - pt + strata(pt) - 1, data = nuclearplants)
  xb2 <- xBalance(pr ~ . - pt + strata(pt) - 1, data = nuclearplants, subset=pt<=1)

  expect_equal(xb1, xb2)
})


test_that("p.adjust.method argument", {
  data(nuclearplants)

  res.none <- xBalance(pr ~ . + strata(pt),
                       data = nuclearplants,
                       report = c("p.value", "chisquare"),
                       p.adjust.method = "none")
  
  # the default argument (holm) should cause the p-values to increase
  res.holm <- xBalance(pr ~ . + strata(pt),
                       data = nuclearplants,
                       report = c("p.value", "chisquare"))

  expect_true(all(res.holm$result[, "p", ] >= res.none$result[, "p", ]))
  expect_true(all(res.holm$overall[, "p.value"] >= res.none$overall[, "p.value"]))
})
