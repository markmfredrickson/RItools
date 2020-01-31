
library('testthat')
context('Inverse assignment weighting')

test_that("Assignment weights properly inferred from treatment assmts by stratum",
          {
tx1 <- c(1,0,1,0,0)
fac1 <- rep(letters[1:2], c(2,3))
fac1 <- as.factor(fac1)
expect_equivalent(ipr(tx1,fac1), c(2,2,3,1.5,1.5))
expect_equivalent(ipr(tx1,as.character(fac1)), c(2,2,3,1.5,1.5))

expect_equivalent(ipr(tx1,fac1,, type="odds vs 1"), c(1,1,1,.5,.5))
}   )


test_that("NA handling for assignment variable",
          {
            ## assignments excluded from assmt probability calcs, receive 0 as ipr wt.
tx2 <- c(1,0,1,0,NA)
fac1 <- rep(letters[1:2], c(2,3))
fac1 <- as.factor(fac1)
expect_equivalent(ipr(tx2,fac1), c(rep(2,4),0))
})

test_that("NA handling for strata variable",
          {
              tx1 <- c(1,0,1,0,0)
              tx2 <- c(1,0,1,0,NA)
              fac2 <- rep(letters[1:2], c(2,3))
              fac2 <- as.factor(fac2)
              fac2[5] <- NA
              expect_equivalent(RItools:::ipr(tx1,fac2), c(rep(2,4),0))
              expect_equivalent(RItools:::ipr(tx2,fac2), c(rep(2,4),0))
          }
          )

test_that("Utility function for specification of reference levels",
          {
              levs <- as.character(0:3)
              expect_error(getRefLev("foo", levs))
              expect_error(getRefLev(c("foo","bar"), levs))
              expect_equivalent(getRefLev("odds vs 1", levs), "1")
              expect_equivalent(getRefLev("odds vs1", levs), "1")
              expect_equivalent(getRefLev("oddsvs1", levs), "1")
              expect_equivalent(getRefLev("oddsvs 1", levs), "1")
              expect_equivalent(getRefLev("odds against 1", levs), "1")

              levs <- c(levs, "01")
              expect_equivalent(getRefLev("odds vs 1", levs), "1")
              expect_equivalent(getRefLev("odds vs 01", levs), "01")
              expect_equivalent(getRefLev("odds against 1", levs), "1")
              expect_equivalent(getRefLev("odds against 01", levs), "01")
          }
          )

test_that("W/ clus arg, wts det'd by cluster- not record-level z * strata table", {
set.seed(201801)
n <- 7L 
dat <- data.frame(y=rnorm(n), x=rnorm(n),
                  s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                  )
dat <- transform(dat, z=as.numeric( (x+rnorm(n))>0 ) )
dat <- transform(dat, inv_pr_wt=ipr(z, s), inv_odds_wt=ipr(z, s, type = "odds vs 1"))
dt <- transform(dat, clus=letters[1L:n])[rep(1L:n, rpois(n, 3)),]
expect_equivalent(dt$inv_pr_wt, with(dt, ipr(z, s, clus)))
expect_equivalent(dt$inv_odds_wt, with(dt, ipr(z, s, clus, type = "odds vs 1")))
})

## To-dos:
## - test dropping of blocks in which not all conditions are represented
