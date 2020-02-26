library("testthat")

context("Stratified Designs")

test_that("Rotated covariates", {
    n <- 30
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)

    x <- matrix(c(x1, x2, x3), ncol = 3)

    s <- factor(c(rep("A", 10), rep("B", n - 10)))
    n1 <- c(A = 3, B  = 5)

    d <- create_stratified_design(s, treated = n1)

    rotx <- rotate_covariates(d, x)

    expect_equal(dim(rotx), c(30, 3))
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
