library('testthat')
library(MASS)

context('Moment Calculation')

## functions for getting empirical moments
source("misc.R")

test_that("Set up code agrees with itself", {

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
        c(1, 1, 1, 1, 1,
          2, 2, 2, 2, 2, 2, 2))
    df$z <- c(1, 0,
              0, 1,
              0, 1, 0,
              0, 1, 0, 1, 1)
    ## end data set up

    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)

    d <- create_stratified_design(df$match, z = df$z)

    emp_t2 <- empirical_t2(x, d@Units, d@Count, d@Treated)
    expect_equal(dim(emp_t2), c(4, 350))

    emp_t2_sums <- colSums(emp_t2)
    expect_equal(mean(emp_t2_sums), 4)

    emp_mal <- empirical_mahalanobis(x, as.matrix(d@Units), d@Count, d@Treated)

    expect_equivalent(emp_t2_sums, emp_mal)

    ## explicit strata level mahalanobis
    rotx <- rotate_covariates(d, x)
    es_S <- empirical_S_by_strata(df$match, c("1" = 2, "2" = 4), rotx)
    indexes <- lapply(es_S, function(x) { 1:dim(x)[2] })
    grd <- as.matrix(do.call(expand.grid, indexes))

    es_S_t <- apply(grd, 1, function(ab) {
        es_S[[1]][, ab[1]] + es_S[[2]][, ab[2]]
    })

    es_S_m <- colSums(es_S_t^2)
    expect_equal(mean(es_S_m), 4)
    expect_equal(var(es_S_m), var(emp_mal))

    expect_equivalent(sort(es_S_m), sort(emp_mal))
})

test_that("Calculating moments of Mahalanobis statistic", {
    library(MASS)

    ### Set up
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

    ### end data set up
    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)
    s <- model.matrix(~ match - 1, data = df)

    ## set up some matrices to indicate strata level stuff
    sn <- as.vector(t(s) %*% rep(1, n))
    sn1 <- as.vector(t(s) %*% df$z)
    sn0 <- sn - sn1

    ## assumes X is in strata sorted order
    zs <- design_to_zs(s, sn, sn1)

    ref_xb_chisquare <- xBalance(z ~ x1 + x2 + x3, data = df, strata = list(a = ~ match), report = 'all')$overall$chisquare

    ms <- empirical_mahalanobis(x, s, sn, sn1)
    expect_equal(mean(ms), 4)

    xbs <- apply(zs, 2, function(z) {
        tmp <- df
        tmp$z <- z
        xb <- xBalance(z ~ x1 + x2 + x3, data = tmp, strata = list(a = ~ match), report = 'all')
        xb$overall$chisquare
    })

    expect_equal(mean(xbs), 4)

    expect_true(all((ms - xbs)^2 <= sqrt(.Machine$double.eps)))


    bts <- apply(zs, 2, function(z) {
        tmp <-df
        tmp$z <- z
        bt <- balanceTest(z ~ x1 + x2 + x3 + strata(match), data = tmp)$overall[1,1]
    })

    expect_equal(mean(bts), 4)

    expect_true(all((ms - bts)^2 <= sqrt(.Machine$double.eps)))
    expect_true(all((xbs - bts)^2 <= sqrt(.Machine$double.eps)))
})

test_that("Direct moment calculations are correct", {
    

### Set up
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

### end data set up
    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)
    s <- model.matrix(~ match - 1, data = df)

    ## set up some matrices to indicate strata level stuff
    sn <- as.vector(t(s) %*% rep(1, n))
    sn1 <- as.vector(t(s) %*% df$z)
    sn0 <- sn - sn1

    ## assumes X is in strata sorted order
    zs <- design_to_zs(s, sn, sn1)

    ### this is largely taken from balanceTest
    fmla <- z ~ x1 + x2 + x3 + strata(match) - 1
    df$`(weights)` <- 1
    design <- makeDesigns(fmla, df)
    aggDesign       <- aggregateDesigns(design)
    aligned <- alignDesignsByStrata(aggDesign)$match

    ## lets double check that we are aligned
    xbar_s <- t(s) %*% x / sn
    xaligned <- x - s %*% xbar_s

    expect_true(all((xaligned - aligned@Covariates[, 1:5])^2 <= sqrt(.Machine$double.eps)))

    ms <- empirical_mahalanobis(x, s, sn, sn1)
    ms_aligned <- empirical_mahalanobis(x, s, sn, sn1)
    expect_true(all(ms == ms_aligned))



})


test_that("Single strata covariance calculations", {
    

    ## Set up
    ## Generate some random data
    ## only one stratum
    set.seed(30303)
    n <- 12
    x1 <- rnorm(n)
    x2 <- x1 + runif(n, -1,  3)
    x3 <- sample(letters[1:3], n, replace = TRUE )
    df <- data.frame(x1, x2, x3)
    df <- df[order(x1), ]
    df$match <- factor((rep("A", 12)))
    df$z <- c(1, 0,
              0, 1,
              0, 1, 0,
              0, 1, 0, 1, 1)
    ## end data set up

    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)

    d <- create_stratified_design(df$match, z = df$z)


    ## compute covariances on natural scale

    t2_cov <- t_squared_covariance(d, x)

    emp_t2 <- empirical_t2(x, d@Units, d@Count, d@Treated)
    emp_mal <- empirical_mahalanobis(x, as.matrix(d@Units), d@Count, d@Treated)

    k <- dim(emp_t2)[2]
    emp_t2_cov <- cov(t(emp_t2)) * (k - 1) / k 

    expect_equal(dim(emp_t2_cov), c(4,4))
    expect_equal(dim(t2_cov), c(4,4))

    expect_equivalent(emp_t2_cov, t2_cov)

    ## Interesting: my hypothesis was wrong. I thought since using the singular
    ## vectors was equivalent to using x for Mahalanobis distance, the
    ## distribution of T^2 would also be the same.
    ## Wrong! The only test that passes compares emp_mal and emp_u_mal
    ## Leaving this in as an interesting little finding.

    ##  xu <- svd(x)$u
    ##  u_cov <- t_squared_covariance(d, xu)
    ##  expect_equivalent(t2_cov, u_cov)

    ##  emp_u <- empirical_t2(xu, d@Units, d@Count, d@Treated)
    ##  k <- dim(emp_u)[2]
    ##  emp_u_cov <- cov(t(emp_u)) * (k - 1) / k 

    ##  emp_u_mal <- empirical_mahalanobis(xu, as.matrix(d@Units), d@Count, d@Treated)
    ##  expect_equivalent(emp_mal, emp_u_mal)

    ##  expect_equivalent(emp_u_cov, emp_t2_cov)
    ##  expect_equivalent(emp_u_cov, u_cov)
})

test_that("Multiple strata covariance calculations", {
    

    ## Set up
    ## Generate some random data
    set.seed(30303)
    n <- 12
    x1 <- rnorm(n)
    x2 <- x1 + runif(n, -1,  3)
    x3 <- sample(letters[1:3], n, replace = TRUE )
    df <- data.frame(x1, x2, x3)
    df <- df[order(x1), ]
    df$match <- factor(c(rep("A", 5),
                         rep("B", 7)))
    df$z <- c(1, 0, 0, 1, 0,
              1, 0, 0, 1, 0, 1, 1)
    ## end data set up

    ## set up the numeric covariance matrix for all the data
    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)

    ## the design using both strata
    d <- create_stratified_design(df$match, z = df$z)
    
    ## now create designs for the A and B strata
    dfa <- df[df$match == "A", ]
    dfa$match <- factor("A")
    da <- create_stratified_design(dfa$match, z = dfa$z)
    
    dfb <- df[df$match == "B", ]
    dfb$match <- factor("B")
    db <- create_stratified_design(dfb$match, z = dfb$z)

    ## generating all possible T^2 and getting cov
    ## helper function first
    covh <- function(tdist) {
        k <- dim(tdist)[2]
        cov(t(tdist)) * (k - 1) / k
    }

    emp_t2 <- empirical_t2(x, d@Units, d@Count, d@Treated)
    emp_t2_cov <- covh(emp_t2)

    emp_t2_a <- empirical_t2(x[df$match == "A", ], da@Units, da@Count, da@Treated)
    emp_t2_a_cov <- covh(emp_t2_a)

    emp_t2_b <- empirical_t2(x[df$match == "B", ], db@Units, db@Count, db@Treated)
    emp_t2_b_cov <- covh(emp_t2_b)


    ## these are the functions in Stratified.R
    t2_cov <- t_squared_covariance(d, x)
    t2_a_cov <- t_squared_covariance(da, x[df$match == "A",])
    t2_b_cov <- t_squared_covariance(db, x[df$match == "B",])

    ## basic checks for number of variables left after rotation
    expect_equal(dim(emp_t2_cov), c(4,4))
    expect_equal(dim(t2_cov), c(4,4))

    ## Now check to make sure everything matches up to the empirical results
    expect_equivalent(emp_t2_a_cov, t2_a_cov)
    expect_equivalent(emp_t2_b_cov, t2_b_cov)
    expect_equivalent(emp_t2_cov, t2_cov)


    ## Breaking down per strata covariance matrices
    ## using eudclidean distance
    emp_strata_S <- empirical_S_by_strata(df$match, c("A" = 2, "B" = 4), x)
    emp_strata_S_cov <- lapply(emp_strata_S, function(st) {
        covh(st)
    })
    emp_strata_S2_cov <- lapply(emp_strata_S, function(st) {
        covh(st^2)
    })

    scm <- strata_t_covariance_matrices(d, x)
    expect_equivalent(scm[1,,], emp_strata_S_cov[[1]])
    expect_equivalent(scm[2,,], emp_strata_S_cov[[2]])

    scm2 <- strata_t2_covariance_matrices(d, x)
    expect_equivalent(scm2[1,,], emp_strata_S2_cov[[1]])
    expect_equivalent(scm2[2,,], emp_strata_S2_cov[[2]])

    ## recompute the covariance of squared euclidean distance
    emp_S2 <- empirical_S2(x, d@Units, d@Count, d@Treated)
    emp_S2_cov <- covh(emp_S2)

    S2_cov <- euclidean_squared_covariance(d, x)

    expect_equivalent(S2_cov, emp_S2_cov)
    expect_equal(sum(S2_cov), sum(emp_S2_cov))
})

test_that("Small strata calculations", {
    

    ## Set up
    ## Generate some random data
    set.seed(30303)
    n <- 12
    x1 <- rnorm(n)
    x2 <- x1 + runif(n, -1,  3)
    x3 <- sample(letters[1:3], n, replace = TRUE )
    df <- data.frame(x1, x2, x3)
    df <- df[order(x1), ]
    df$match <- factor(
        c("A", "A",
          "B", "B",
          "C", "C", "C",
          "D", "D", "D",
          "E", "E"))
    df$z <- c(1, 0,
              0, 1,
              0, 1, 0,
              0, 1, 1,
              1, 0)
    ## end data set up

    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)
    d <- create_stratified_design(df$match, z = df$z)

    ## generating all possible T^2 and getting cov
    ## helper function first
    covh <- function(tdist) {
        k <- dim(tdist)[2]
        cov(t(tdist)) * (k - 1) / k
    }

    emp_t2 <- empirical_t2(x, d@Units, d@Count, d@Treated)
    emp_t2_cov <- covh(emp_t2)
    t2_cov <- t_squared_covariance(d, x)
    expect_equivalent(emp_t2_cov, t2_cov)

})
