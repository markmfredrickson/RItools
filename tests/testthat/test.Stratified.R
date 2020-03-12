library("testthat")

context("Stratified Designs")

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

    ## a version of this function also appears in test.moments.R
    empirical_cov <- function(x, s, sn, sn1) {
        sn0 <- sn - sn1

        zs <- design_to_zs(s, sn, sn1)
        k <- dim(zs)[2]

        ## strata level probabilities
        p <- sn1 / sn # pi

        ZtoJ <- function(z) { (z - s %*% p) / (s %*% p) }

        ## Compute the Mahalanobis distance
        xj <- function(z) {
            j <- ZtoJ(z)
            t(x) %*% j
        }

        xjs <- apply(zs, 2, xj)

        cov_xj <- (k - 1) / k * cov(t(xjs))

        return(cov_xj)
    }


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

