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

    expect_equal(as.matrix(sdz@Units), sm)
    expect_equal(as.matrix(sdn@Units), sm)

    expect_equal(sdz@Count, c(2, 3))
    expect_equal(sdn@Count, c(2, 3))

    expect_equal(sdn@Treated, c(1, 2))
    expect_equal(sdn@Treated, c(1, 2))
})
