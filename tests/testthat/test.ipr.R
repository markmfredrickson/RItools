
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


## To-dos:
## - test cluster arg
## - test dropping of blocks in which not all conditions are represented
