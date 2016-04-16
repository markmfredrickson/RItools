################################################################################
# Tests for xBalance function
################################################################################
## if working interactively in inst/tests you'll need
## library(RItools, lib.loc = '../../.local')
library("testthat")

context("xBalance Functions")

test_that("xBal univariate descriptive means agree w/ reference calculations",{
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
     expect_equivalent(xb1$results["x1", "adj.diff", "s"], coef(lm2a)["z"])

    ## a little more explicitly:
    d <- split(dat, dat$s)
    mndiffs <- sapply(d, function(Data) {with(Data, mean(x1[z==1]) - mean(x1[z==0]))})
    cmndiff <- weighted.mean(mndiffs, w = sapply(d, function(Data) sum(Data$z==1)))
    expect_equivalent(xb1$results["x1", "adj.diff", "s"], cmndiff)

})

test_that("xBal univariate inferentials, incl. agreement w/ Rao score test for cond'l logistic regr",{
    library(survival)
    set.seed(20160406)
    n <- 7 # increase at your peril -- clogit gets slow quickly as stratum size increases
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

test_that("Alternate formats for stratum.weights argument", {
    set.seed(20160406)
    n <- 7 # increase at your peril -- clogit gets slow quickly as stratum size increases
    dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                      s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                      )
    dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )

    xb1 <- xBalance(z~x1+strata(s)-1, data=dat, report="all")

    hwts <- with(dat, colSums(table(z, s)^-1)^-1 ) # 2*harmonic means of (n_{tb}, n_{cb}), not normalized
    xb1a <- xBalance(z~x1+strata(s)-1, data=dat, stratum.weights=hwts, report="all")
    expect_identical(xb1, xb1a)

    xb2 <- xBalance(z~x1+strata(s), data=dat, report="all")
    xb2a <- xBalance(z~x1+strata(s), data=dat, stratum.weights=list(Unstrat=c("1"=1), s=hwts), report="all")
    expect_identical(xb2, xb2a)
    xb2b <- xBalance(z~x1+strata(s), data=dat, stratum.weights=list(Unstrat=1, s=hwts), report="all")
    expect_identical(xb2, xb2b)
    xb2c <- xBalance(z~x1+strata(s), data=dat,
                     stratum.weights=list(Unstrat="cheese!", #shouldn't matter in 1-stratum case
                                                   s=hwts), report="all")
    expect_identical(xb2, xb2c)
    xb2d <- xBalance(z~x1+strata(s), data=dat,
                     stratum.weights=list(Unstrat=NULL, s=hwts), report="all")
    expect_identical(xb2, xb2d)
    xb2e <- xBalance(z~x1+strata(s), data=dat,
                     stratum.weights=list(s=hwts), report="all")
    expect_identical(xb2, xb2e)

} )
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

  ## Developer note: to strip out entries corresponding to intercept -- which has var 0,
  ## except when there's variation in element weights and/or cluster sizes --
  ## have to filter out rows and cols named "(element weight)", separately for each
  ## entry in list tcov.  (Recording while updating test that follows, `c(4,4)` --> `c(5,5)`)
  expect_equal(dim(tcov[[1]]), c(5,5))

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
  res6 <- xBalance(pr ~ cost, data=nuclearplants, post.alignment.transform = rank,
                   report="all", p.adjust.method='none')

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

  xb1 <- xBalance(pr ~ . - pt + strata(pt) - 1, data = nuclearplants)
  xb2 <- xBalance(pr ~ . - pt + strata(pt) - 1, data = nuclearplants, subset=pt<=1)
  expect_equal(xb1, xb2)

  n2 <- nuclearplants
  n2 <- rbind(n2, n2[1,])
  n2[nrow(nuclearplants)+1, "pt"] <- 2

  xb3 <- xBalance(pr ~ . - pt + strata(pt) - 1, data = n2, subset=pt<=1)
  expect_equal(xb1, xb3)
})

test_that("Observations not meeting subset condition are retained although downweighted to 0",{

    data(nuclearplants)
    ## first, check assumptions about offsets that are made within the code
    mf0 <- model.frame(cost~date + offset(date<68), data=nuclearplants, offset=(cap>1000))
    expect_equal(sum(names(mf0)=='(offset)'), 1L)
    expect_equivalent(mf0$'(offset)', nuclearplants$cap>1000)
    
    n2 <- nuclearplants
    nuclearplants$pt <- factor(nuclearplants$pt)
    n2 <- rbind(n2, n2[1,])
    n2[nrow(nuclearplants)+1, "pt"] <- 2
    n2$pt <- factor(n2$pt)

    ## this indirect test relies on xBal's dropping unused factor levels
    xb1 <- xBalance(pr ~ ., data = nuclearplants)
    xb2 <- xBalance(pr ~ ., data = n2, subset=pt!='2')
    ## confirm that we still see the '2' level, even if it receives no weight
    expect_match(dimnames(xb2$results)[[1]], "pt2", all=FALSE)
    expect_equivalent(xb1$results[,'std.diff',], #only the descriptives should be the same for
                      xb2$results[ dimnames(xb2$results)[[1]]!="pt2" ,'std.diff',]) #these two
    expect_true(is.na(xb2$results[ "pt2" ,'std.diff',]) | # presently this is NA, but
                    xb2$results[ "pt2" ,'std.diff',]==0) # it might ideally be a 0
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
  T_or_NA <- function(vec) {ans <- as.logical(vec) ; ans[is.na(ans)] <- TRUE; ans}
  
  expect_true(all(T_or_NA(res.holm$result[, "p", ] >= res.none$result[, "p", ])))
  expect_true(all(T_or_NA(res.holm$overall[, "p.value"] >= res.none$overall[, "p.value"])))

### with just one covar, holm should do the same as none

    res1.none <- xBalance(pr ~ cost + strata(pt),
                       data = nuclearplants,
                       report = c("p.value", "chisquare"),
                       p.adjust.method = "none")
  
  res1.holm <- xBalance(pr ~ cost + strata(pt),
                       data = nuclearplants,
                       report = c("p.value", "chisquare"))

  expect_equal(res1.holm$result[, "p", ], res1.none$result[, "p", ])
  expect_equal(res1.holm$overall[, "p.value"], res1.none$overall[, "p.value"])

})

## To do: adapt the below to test print.xbal instead of lower level functions
##test_that("printing of NA comparisons is optional",
replicate(0,
{
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

  design.flags   <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d)
  design.noFlags <- RItools:::makeDesigns(z ~ x + f + strata(s), data = d, include.NA.flags = FALSE)

  expect_equal(dim(design.flags@Covariates)[2], 5)
  expect_equal(dim(design.noFlags@Covariates)[2], 4)
})
