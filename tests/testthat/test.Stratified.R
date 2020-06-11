library("testthat")

context("Stratified Designs")

## A copy of this also exists in test.Stratified.R
## TODO: Move them to own utils.R file

## @param s A n x s stratum membership matrix
## @param sn A s-vector of units per stratum
## @param sn1 A s-vector of treated units per stratum
## @return A matrix of treatment assignments at the unit level
design_to_zs <- function(s, sn, sn1) {

    ## generate all possible randomizations
    zs_by_strata <- mapply(sn, sn1, FUN = function(n, n1) {
        apply(combn(n, n1), 2, function(idx) {
            zz <- numeric(n)
            zz[idx] <- 1
            return(zz)
        })
    }, SIMPLIFY = FALSE)

    zj <- sapply(zs_by_strata, ncol)
    idxes <- as.matrix(do.call(expand.grid, lapply(zj, function(k) { 1:k })))
    zs <- apply(idxes, 1, function(idx) {
        unlist(mapply(idx, zs_by_strata, FUN = function(a, b) { b[, a] }))})

    return(zs)
}

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

        ## j = Z - E(Z)
        ZtoJ <- function(z) { (z - s %*% p) }

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
