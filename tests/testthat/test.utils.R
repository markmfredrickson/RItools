################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("Utility Functions")

test_that("xbal formula method", {
  # create a quick balanceTest result
  df <- data.frame(Z = rep(c(1,0), 10),
                   X = rnorm(20))

  xb <- balanceTest(Z ~ X, data = df)

  # we expect to have a no intercept formula
  expect_equal(Z ~ X, as.formula(xb))
})


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

              slm.n1 <- slm_fit_csr(nullfac.csr, quickY)

              expect_equal(lm.n$fitted, as.vector(slm.n1$fitted))
              expect_equal(lm.n$residuals, as.vector(slm.n1$residuals))
          }
          )
test_that("sparse design strat mean calculator returns 0 for strata w/o non-null wts",
          {
              expect_true(require("SparseM"))
              fac <- factor(rep(c("a", "b"), each=2))
              wts  <- c(1:2, 0, 0)
              S_ <- SparseMMFromFactor(fac)
              # .lm.fit() handles singular fits without warnings
              #expect_warning(theslmfit  <- slm.wfit.csr(S_, matrix(1:4), weights=wts),
              #               "singularity problem")
              theslmfit  <- slm.wfit.csr(S_, matrix(1:4), weights=wts)
              expect_equivalent(theslmfit$fitted,
                                c(rep((1*1+2*2)/(1+2), 2), 0, 0))
          }
)

test_that("gramian_reduction works as expected",
          {
            expect_true(require("SparseM"))
            #let tl stand in for xprimex
            tl <- diag(1, 5)
            diag(tl)[c(1,2,4)] <- 0
            zeroes <- diag(tl) == 0
            not.zeroes <- !zeroes
            stl <- tl[, not.zeroes]
            sparse_matrix <- gramian_reduction(zeroes)
            
            expect_identical(stl, as.matrix(sparse_matrix))
            
            tl <- diag(1, 5)
            diag(tl)[c(2,4)] <- 0
            zeroes <- diag(tl) == 0
            not.zeroes <- !zeroes
            stl <- tl[, not.zeroes]
            sparse_matrix <- gramian_reduction(zeroes)
            
            expect_identical(stl, as.matrix(sparse_matrix))
            
            
            tl <- diag(1, 5)
            diag(tl)[c(2)] <- 0
            zeroes <- diag(tl) == 0
            not.zeroes <- !zeroes
            stl <- tl[, not.zeroes]
            sparse_matrix <- gramian_reduction(zeroes)
            
            expect_identical(stl, as.matrix(sparse_matrix))
            
            
            tl <- diag(1, 5)
            diag(tl)[c(5)] <- 0
            zeroes <- diag(tl) == 0
            not.zeroes <- !zeroes
            stl <- tl[, not.zeroes]
            sparse_matrix <- gramian_reduction(zeroes)
            
            expect_identical(stl, as.matrix(sparse_matrix))
            
            tl <- diag(1, 5)
            diag(tl)[c(1, 5)] <- 0
            zeroes <- diag(tl) == 0
            not.zeroes <- !zeroes
            stl <- tl[, not.zeroes]
            sparse_matrix <- gramian_reduction(zeroes)
            
            expect_identical(stl, as.matrix(sparse_matrix))
            
            tl <- diag(1, 2)
            diag(tl)[c(1)] <- 0
            zeroes <- diag(tl) == 0
            not.zeroes <- !zeroes
            stl <- tl[, not.zeroes, drop = FALSE]
            sparse_matrix <- gramian_reduction(zeroes)
            
            expect_identical(stl, as.matrix(sparse_matrix))
            
            tl <- diag(1, 2)
            diag(tl)[c(1, 2)] <- 0
            zeroes <- diag(tl) == 0
            
            expect_error(gramian_reduction(zeroes)) 
            
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

test_that("Select a subset of xbal results (for printing, etc)", {

### create a quick balanceTest result
  set.seed(20121129)
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = cut(rnorm(n), breaks = 3),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   U = rnorm(n),
                   S = as.factor(rep(c("A", "B"), each = n/2)))


  # the data are will be grouped as (X, XY) and (all levels W by all levels
  # K), Y is not in an explicit group
  xb <- balanceTest(Z ~ X * Y + W * K + strata(S) + strata(cut(U, 3)),
                 data = df)

  # strata based subsetion is the easiest as it common across groups and
  # variables
  xb.noneU <- subset(xb, strata = c("--", "cut(U, 3)"))

  expect_equivalent(dim(xb$results)[1:2], dim(xb.noneU$results)[1:2])
  expect_equivalent(dim(xb$results)[3]-1, dim(xb.noneU$results)[3])
  expect_equal(attr(xb$results, "originals"),
               attr(xb.noneU$results, "originals"))

  # stat and test subsetion only limits the results and overall tables
  # repectively

  xb.zp <- subset(xb, stats = c("z", "p"))
  expect_equal(attr(xb$results, "originals"),
               attr(xb.zp$results, "originals"))
  expect_equal(xb$overall, xb.zp$overall)
  expect_equivalent(dim(xb.zp$results)[c(1,3)], dim(xb$results)[c(1,3)])
  expect_equivalent(dim(xb.zp$results)[2], 2)

  xb.chip <- subset(xb, tests = c("chisquare", "p.value"))
  expect_equal(xb.chip$results, xb$results)
### don't need next one if prev expectation checks out
###  expect_equal(attr(xb$results, "originals"),
###               attr(xb.chip$results, "originals"))
  expect_equivalent(dim(xb.chip$overall), c(3, 2))

  # for vars and groups, we might want the arguments to interact more directly
  # e.g. subseting vars only returns groups that include the vars
  # e.g. subseting groups only returns the variables contained in those groups

  # for vars, only limit the results table (for now)
  # use exact variable names -- possibly use regexes later
  xb.XY <- subset(xb, vars = c("X", "Y"))

  expect_equivalent(dim(xb.XY$results), c(2, 7, 3))
  expect_equal(length(attr(xb.XY$results, "originals")), 2)


})

test_that("Formatting w/ appropriate sigfigs values",{

  ## same start point as last test
  set.seed(20121129)
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = cut(rnorm(n), breaks = 3),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   U = rnorm(n) )
  xb1 <- balanceTest(Z ~ X * Y + W * K + strata(cut(U, 3)),
                 data = df)
  pt1 <- print(xb1, digits=2, printme=FALSE)
  xb2 <- balanceTest(Z ~ X * Y + W * K + strata(cut(U, 3)),
                 data = df)
  ## Confirm that the by-row rounder rounds to requested sigfigs, potentially
  ## w/ more rounding for Control and Treatment columns:
  ca1 <- original_units_var_formatter(xb2$results[,c("Control", "Treatment", "adj.diff"),, drop=F], 3)
  expect_equal(ca1[1,,1], c("Control"="-0.0474",  "Treatment"="-0.0425",  "adj.diff"=" 0.00495"))

  pt2 <- print(xb2, digits=2, printme=FALSE)


  ## Brittle test: with settings as above, I'm pretty sure 2 sigfigs
  ## should translate to the thousandths place for the z-stats.
  expect_equivalent(round(as.numeric(pt1[[1]][,4]), 4), as.numeric(pt1[[1]][,4]))
  ## This similarly brittle pair of tests confirms that rounding is being done
  ## separately for different variables, at least as it affects columns
  ## Treatment and Control
  expect_equivalent(round(as.numeric(pt1[[1]][3,1:2]), 2), as.numeric(pt1[[1]][3,1:2]))
  expect_gt(abs( round(as.numeric(pt1[[1]][1,1]), 2) - as.numeric(pt1[[1]][1,1]) ),0)
  expect_equivalent(round(as.numeric(pt2[[1]][3,1:2]), 2), as.numeric(pt2[[1]][3,1:2]))
  expect_gt(abs( round(as.numeric(pt2[[1]][1,1]), 2) - as.numeric(pt2[[1]][1,1]) ),0)

})

test_that("Var names in print.xbal(,horiz=F)",{
  set.seed(20121129)
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = cut(rnorm(n), breaks = 3),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   U = rnorm(n),
                   S = as.factor(rep(c("A", "B"), each = n/2)))


  # the data are will be grouped as (X, XY) and (all levels W by all levels
  # K), Y is not in an explicit group
  xb <- balanceTest(Z ~ X * Y + W * K + strata(S) + strata(cut(U, 3)),
                 data = df)

  pt2 <- print(xb, horiz=F, printme=F)
  expect_true('std.diff' %in% colnames(pt2$vartable[['S']])) # not 'std.diff.S'
})

test_that("Matrix pseudoinversion also returns rank", {
  m1 <- diag(c(1,1,0))
  a1 <- RItools:::XtX_pseudoinv_sqrt(m1)
  expect_equal(attr(a1, "r"), 2)


  x <- rnorm(100)
  xx <- matrix(c(x, x, x), ncol = 3)
  xtx <- t(xx) %*% xx
  a2 <- RItools:::XtX_pseudoinv_sqrt(xtx)
  expect_equal(attr(a2, "r"), 1)



})

