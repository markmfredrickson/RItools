library("testthat")

context("Stratified Designs")

source("misc.R")

test_that("Rotated covariates", {
    n <- 12 
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)

    x <- matrix(c(x1, x2, x3), ncol = 3)

    s <- factor(c(rep("A", n / 2), rep("B", n / 2)))
    n1 <- c(A = n/4, B  = n/4 + 1)

    d <- create_stratified_design(s, treated = n1)

    rotx <- rotate_covariates(d, x)

    expect_equal(dim(rotx), c(12, 3))

    emp_cov <- empirical_cov(x, d@Units, d@Count, d@Treated)
    emp_rot <- XtX_pseudoinv_sqrt(emp_cov, TRUE)
    expect_equivalent(rotx@Rotation, emp_rot)
    expect_equivalent(x %*% emp_rot, as.matrix(rotx@.Data))
})

test_that("Lower rank X", {

    n <- 30
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- x1 + x2

    x <- matrix(c(x1, x2, x3), ncol = 3)

    s <- factor(c(rep("A", 10), rep("B", n - 10)))
    n1 <- c(A = 3, B  = 5)

    d <- create_stratified_design(s, treated = n1)

    rotx <- rotate_covariates(d, x)

    expect_equal(dim(rotx), c(30, 2))

    xc <- sample(c("A", "B", "C"), n, replace = TRUE)
    df <- data.frame(x1, x2, xc)
    xx <- model.matrix(~ x1 + x2 + xc - 1, df)
    
    rotxx <- rotate_covariates(d, xx)
    expect_equal(dim(rotx), c(30, 2))
})


test_that("Basic stratified designs using either counts or Z", {
    s <- factor(c("A", "B", "A", "B", "B"))
    z <- c(T, T, F, F, T)
    sn1 <- c(B = 2, A = 1)

    sdz <- create_stratified_design(s, z = z)
    sdn <- create_stratified_design(s, treated = sn1)

    sm <- matrix(c(1, 0, 1, 0, 0,
                   0, 1, 0, 1, 1), ncol = 2)

    expect_equivalent(as.matrix(sdz@Units), sm)
    expect_equivalent(as.matrix(sdn@Units), sm)

    expect_equal(sdz@Count, c(2, 3))
    expect_equal(sdn@Count, c(2, 3))

    expect_equal(sdn@Treated, c(1, 2))
    expect_equal(sdn@Treated, c(1, 2))
})

test_that("balanceTest method for Stratified objects", {
  ## Set up
  ## Generate some random data
  set.seed(30303)
  n <- 12
  x1 <- rnorm(n)
  x2 <- x1 + runif(n, -1,  3)
  x3 <- sample(letters[1:3], n, replace = TRUE )
  df <- data.frame(x1, x2, x3)
  df <- df[order(x1), ]
  df$match <- as.factor(
      c(1, 1,
        2, 2,
        3, 3, 3,
        4, 4, 4, 4, 4))
  df$z <- c(1, 0,
            0, 1,
            0, 1, 0,
            0, 1, 0, 1, 1)

  ## end data set up
  x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)

  btf <- balanceTest(z ~ x1 + x2 + x3 + strata(match), data = df)

  ## making the stratified object directly
  strat <- create_stratified_design(strata = df$match, z = df$z)

  bts <- balanceTest(strat, data = x, z = df$z)
  expect_is(bts, "list")
  expect_is(bts[[1]], "MahalanobisDistance")
  expect_is(bts[[1]], "BalanceTest")

  expect_equivalent(bts, btf[1]) ## just pull the first item (as a list) from the formula version

  ## explicit .method names to make sure we are using the right method dispatch in the main function.
  rot <- rotate_covariates.StratifiedDesign(strat, x)
  dst <- mahalanobis_distance.DesignRotatedCovariates(rot, df$z)
  expect_equal(bts[[1]]@.Data, dst@.Data)

  pv <- pvalue.SecondOrderChisquareApproximation(dst@Distribution, dst@.Data)
  expect_equal(bts[[1]]@Pvalue, pv)


  ## just for completeness
  unstrat <- create_stratified_design(strata = as.factor(rep(1, n)), z = df$z)
  btu <- balanceTest(unstrat, data = x, z = df$z)
  expect_equal(btu[[1]], btf[[2]])

})

test_that("Pairwise products and strata means", {
  ## tests for internal functions 
  ## Stratified.R/pairwise_products
  ## Stratified.R/strata_pairwise_means
  
  # n = 3, k = 2
  a <- matrix(0:5, ncol = 2)
  aa <- RItools:::pairwise_products(a)
  
  # expect dim of (n, k, k)
  expect_equal(dim(aa), c(3, 2, 2))
  
  expected_aa <- array(c(0, 1, 4, 0, 4, 10,
                         0, 4, 10, 9, 16, 25), dim = c(3, 2, 2))
  
  expect_equal(aa, expected_aa)
                       
})

test_that("Switching to Matrix class", {
  
  ## Tests to verify that our matrix operations survive a change to the Matrix class
  v <- factor(c("c", "a", "a", "c", "b", "a"))
  xx <- MatrixFromFactor(v)
  
  expect_equal(dim(xx), c(6, 3))
  
  u <- matrix(c(0, 0, 1,
                1, 0, 0,
                1, 0, 0,
                0, 0, 1,
                0, 1, 0,
                1, 0, 0), byrow = TRUE, ncol = 3)
  
  expect_equivalent(as.matrix(xx), u)
})
