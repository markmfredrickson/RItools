
library('testthat')
context('Inverse assignment weighting')

test_that("Assignment weights properly inferred from treatment assmts by stratum",
          {
tx1 <- c(1,0,1,0,0)
fac1 <- rep(letters[1:2], c(2,3))
fac1 <- as.factor(fac1)
expect_equivalent(ipwts_blockRCT(tx1,fac1), c(2,2,3,1.5,1.5))
expect_equivalent(ipwts_blockRCT(tx1,as.character(fac1)), c(2,2,3,1.5,1.5))
}   )


test_that("NA handling for assignment variable",
          {
            ## assignments excluded from assmt probability calcs, receive 0 as ipw.
tx2 <- c(1,0,1,0,NA)
fac1 <- rep(letters[1:2], c(2,3))
fac1 <- as.factor(fac1)
expect_equivalent(ipwts_blockRCT(tx2,fac1), c(rep(2,4),0))
          })

## To-dos:
## - test cluster arg
## - test dropping of blocks in which not all conditions are represented
