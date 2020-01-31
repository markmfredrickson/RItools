################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("Utility Functions")



test_that("fitter for sparse designs handles intercept only design",
          {
              expect_true(require("SparseM"))
              nullfac <- factor(rep("a", 4))
              nullfac.csr <- SparseMMFromFactor(nullfac)
              quickY <- 1:4
              quickY.mat <- matrix(quickY)
              quickY.csr <- as(quickY.mat, "matrix.csr")
              expect_equivalent(SparseM::slm.fit.csr(nullfac.csr, quickY), # i.e., no problem
                                SparseM::slm.fit.csr(nullfac.csr, quickY.mat) # w/ these y classes
                                )
              expect_error(SparseM::slm.fit.csr(nullfac.csr, quickY.csr)) # when this ceases,
                                        # we get to revert to SparseM:slm.fit.csr
              lm.n <- lm.fit(matrix(1,4,1), quickY)

              slm.n1 <- slm.fit.csr.fixed(nullfac.csr, quickY)

              expect_equal(lm.n$fitted, as.vector(slm.n1$fitted))
              expect_equal(lm.n$residuals, as.vector(slm.n1$residuals))
          }
          )

test_that("Residuals from weighted regressions w/ sparse designs",
          {
              nullfac <- factor(rep("a", 4))
              nullfac.csr <- SparseMMFromFactor(nullfac)

              quickfac <- factor(rep(letters[1:2], each=2))
              quickfac.csr <- SparseMMFromFactor(quickfac)
              quickY <- 1:4

              ## First, document the problem.
              slm.u <- SparseM:::slm.fit(quickfac.csr, quickY)
              slm.w <- SparseM:::slm.wfit(quickfac.csr, quickY,
                                         weights=rep(1:2, each=2)
                                         )
              expect_equal(coef(slm.u), # this 
                           coef(slm.w)) # is OK...

              ## but since SparseM:::slm.wfit() doesn't compensate for weighting of design matrix:
              expect_false(all(abs(slm.u$resid - slm.w$resid) < .Machine$double.eps^.5))
              expect_false(all(abs(slm.u$fitted - slm.w$fitted) < .Machine$double.eps^.5))

              ## Now for the fix
              ## First, example from `lm.wfit`:
              set.seed(129)
              
              n <- 7 ; p <- 2
              X <- matrix(rnorm(n * p), n, p) # no intercept!
              y <- rnorm(n)
              w <- rnorm(n)^2
     
              lmw <- lm.wfit(x = X, y = y, w = w)

              X.csr <- as(X, "matrix.csr")

              slmw <- slm.wfit.csr(x=X.csr, y=y, weights=w)

              expect_equal(coef(lmw), as.vector(coef(slmw)), check.attributes=FALSE) # OK here...

              expect_equal(lmw$fitted, as.vector(slmw$fitted)) # OK here 
              expect_equal(lmw$residuals, as.vector(slmw$residuals)) # also

          })


context("Printing of results")



