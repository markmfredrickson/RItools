library('testthat')
context('Moment Calculation')

test_that("Calculating moments of Mahalanobis statistic", {
    library(MASS)

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

    x <- model.matrix(~ x1 + x2 + x3 - 1, data = df)
    s <- model.matrix(~ match - 1, data = df)

    ## set up some matrices to indicate strata level stuff
    sn <- as.vector(t(s) %*% rep(1, n))
    sn1 <- as.vector(t(s) %*% df$z)
    sn0 <- sn - sn1

    ## strata level probabilities
    p <- sn1 / sn # pi
    one_p <- sn0 / sn # 1 - pi
    p_one_p <- p * one_p # pi (1 - pi)

    ## The E(JJ')
    Gamma <- s %*% diag((p * (sn1 - 1) / (sn - 1) - p^2) / (p^2 * one_p^2))  %*% t(s)
    diag(Gamma) <- 1 / as.vector(s %*% p_one_p)

    ## now combine Gamma with x get the central matrix of the treatment assignment sandwich
    x_svd <- svd(x)
    u <- x_svd$u[, x_svd$d >= sqrt(.Machine$double.eps)]

    G_inv <- ginv(Gamma)

    uuG_uu <- u %*% t(u) %*% G_inv %*% u %*% t(u)

    ## Compute the Mahalanobis distance
    mahal <- function(z) {
        j <- (z - s %*% p) / (s %*% p_one_p)
        (t(j) %*% uuG_uu %*% j)[1,1]
    }

    mahal2 <- function(z) {
        j <- (z - s %*% p) / (s %*% p_one_p)
        (t(j) %*% x %*% ginv(t(x) %*% Gamma %*% x) %*% t(x) %*% j)[1,1]
    }

    ref_xb_chisquare <- xBalance(z ~ x1 + x2 + x3, data = df, strata = list(a = ~ match), report = 'all')$overall$chisquare
    expect_equal(mahal(df$z),
                 ref_xb_chisquare)

    expect_lte((mahal2(df$z) - ref_xb_chisquare)^2, sqrt(.Machine$double.eps))

    ## generate all possible randomizations
    zs_by_strata <- mapply(sn, sn1, FUN = function(n, n1) {
       apply(combn(n, n1), 2, function(idx) {
            zz <- numeric(n)
            zz[idx] <- 1
            return(zz)
        })
    })

    js_by_strata <- mapply(1:nlevels(df$match), zs_by_strata, FUN = function(id, zs) {
        js <- (zs - p[id]) / p_one_p[id]
    })

    nj <- sapply(js_by_strata, ncol)
    idxes <- as.matrix(do.call(expand.grid, lapply(nj, function(k) { 1:k })))

    ms <- apply(idxes, 1, function(idx) {
        ## BEWARE makes strong assumption uuG_uu is in strata sorted order
        j <- unlist(mapply(idx, js_by_strata, FUN = function(a, b) { b[, a]}))
        t(j) %*% uuG_uu %*% j
    })

    expect_equal(mean(ms), 4)

    xbs <- apply(idxes, 1, function(idx) {
        tmp <- df
        ## Same assumption as above
        tmp$z <- unlist(mapply(idx, zs_by_strata, FUN = function(a, b) { b[, a]}))
        xb <- xBalance(z ~ x1 + x2 + x3, data = tmp, strata = list(a = ~ match), report = 'all')
        xb$overall$chisquare
    })

    expect_equal(mean(xbs), 4)

    ms2 <- apply(idxes, 1, function(idx) {
        z <- unlist(mapply(idx, zs_by_strata, FUN = function(a, b) { b[, a]}))
        mahal2(z)
    })

    expect_equal(mean(ms2), 4)

    expect_true(all((ms2 - xbs)^2 <= sqrt(.Machine$double.eps)))
})
