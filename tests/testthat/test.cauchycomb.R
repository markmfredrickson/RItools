################################################################################
## Tests for the Cauchy combination (ACAT) omnibus balance statistic.
##
## Motivation (the substantive point these tests defend):
##
## The HB08 d^2 omnibus inverts the p-by-p permutation covariance of the
## adjusted differences.  When the number of covariates p reaches the
## within-block residual degrees of freedom (N - B), that covariance is
## rank-deficient, the statistic freezes at its rank, and the chi-square p-value
## freezes near 0.46 -- the covariates stop mattering.  See the high-dimensional
## section of vignettes/resolution-profile-memo.qmd.
##
## ACAT (Liu and Xie 2020, JASA 115:393-402) sidesteps this: it COMBINES the
## per-covariate p-values instead of inverting their covariance, so it has no
## rank ceiling and stays valid when p exceeds N.  The combination is
##
##     T = sum_k w_k * tan((0.5 - p_k) * pi),   p_ACAT = 0.5 - atan(T)/pi,
##
## the survival function of a standard Cauchy.  Validity rests on each p_k being
## uniform under the null and on the Cauchy's stability under averaging (which is
## what makes it robust to dependence among the covariates).
##
## These tests encode the statistical principles, not just "runs without error":
## the published formula (idempotence), null calibration (uniformity), and the
## design point (sensitivity to a single strong signal; no high-dimensional
## collapse).  Helper is the internal RItools:::acat_pvalue.
################################################################################
library("testthat")
context("cauchy combination (ACAT) omnibus")

## --- 1.  The published formula, checked by idempotence --------------------
##
## Combining d identical p-values must return that p-value exactly.  With all
## p_k = q and equal weights, T = tan((0.5 - q)*pi), and since
## atan(tan(theta)) = theta for theta = (0.5 - q)*pi in (-pi/2, pi/2) whenever
## q is in (0, 1), we get 0.5 - atan(T)/pi = q.  This pins down BOTH the
## transform and the final survival step; the buggy branch version (which halves
## p and then doubles a min) fails it.

test_that("acat_pvalue reproduces an identical p-value (idempotence)", {
  expect_equal(RItools:::acat_pvalue(rep(0.30, 5)), 0.30, tolerance = 1e-8)
  expect_equal(RItools:::acat_pvalue(rep(0.73, 3)), 0.73, tolerance = 1e-8)
  expect_equal(RItools:::acat_pvalue(rep(0.05, 8)), 0.05, tolerance = 1e-8)
})

## --- 2.  Null calibration: a combined p-value must be uniform --------------
##
## Under the global null every per-covariate p_k is Uniform(0,1).  A valid
## combination is then itself Uniform(0,1): mean 0.5 and rejection rates equal to
## their nominal levels.  This is the property that justifies using the number as
## a p-value at all, so it is the most important test here.

test_that("acat_pvalue is calibrated under the global null", {
  skip_on_cran()                       # Monte Carlo; deterministic via the seed
  set.seed(2024)
  B <- 4000; K <- 6
  null_p <- replicate(B, RItools:::acat_pvalue(runif(K)))
  expect_equal(mean(null_p), 0.5, tolerance = 0.02)
  expect_equal(mean(null_p < 0.05), 0.05, tolerance = 0.015)
  expect_equal(mean(null_p < 0.10), 0.10, tolerance = 0.02)
  expect_gt(suppressWarnings(ks.test(null_p, "punif")$p.value), 0.01)
})

## --- 3.  Sensitivity: one strong signal drives the combination ------------
##
## The design point.  ACAT is meant to flag imbalance carried by even a single
## covariate, so a tiny p_k among otherwise-large ones must pull the combined
## p-value small.  (This is the behavior the rank-deficient d^2 cannot deliver in
## high dimension, where it is pinned near 0.46.)

test_that("a single tiny p-value drives the combination toward significance", {
  expect_lt(RItools:::acat_pvalue(c(1e-8, 0.6, 0.7, 0.8, 0.9)), 0.05)
  ## and conversely, uniformly-large p-values do not manufacture significance
  expect_gt(RItools:::acat_pvalue(c(0.6, 0.7, 0.8, 0.9)), 0.30)
})

## --- 4.  Robustness of the implementation ---------------------------------
##
## Boundary p-values (0 or 1) send tan() to +/- Inf; the function must guard
## them rather than returning NaN.  NA p-values (e.g. a covariate with no
## within-block variation) must be dropped, not poison the whole combination.
## The result must always be a valid probability.

test_that("acat_pvalue handles boundaries, NAs, and stays in (0,1)", {
  expect_true(is.finite(RItools:::acat_pvalue(c(0, 0.5, 1))))
  expect_equal(RItools:::acat_pvalue(c(0.3, NA, 0.3)), 0.30, tolerance = 1e-8)
  expect_true(is.na(RItools:::acat_pvalue(c(NA, NA))))
  set.seed(7)
  for (i in 1:50) {
    val <- RItools:::acat_pvalue(runif(sample(2:8, 1)))
    expect_true(val > 0 && val < 1)
  }
  ## default (equal) weights agree with explicitly equal weights
  expect_equal(RItools:::acat_pvalue(c(0.1, 0.5, 0.9)),
               RItools:::acat_pvalue(c(0.1, 0.5, 0.9), weights = rep(1, 3)),
               tolerance = 1e-10)
})

## --- 5.  Integration: opt-in, and off by default --------------------------
##
## The feature mirrors sigma_x_test: it must be requested.  Off by default the
## $overall table is byte-for-byte what it was, which protects the legacy
## print snapshots (tests/*.Rout.save).

test_that("cauchy.combination is opt-in and leaves default output unchanged", {
  set.seed(1); n <- 40
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                  z = rep(c(1, 0), n / 2), s = factor(rep(1:4, each = n / 4)))
  bt_default <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = d)
  bt_acat    <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = d,
                            cauchy.combination = TRUE)
  expect_false("cauchy_comb_p" %in% colnames(bt_default$overall))
  expect_true("cauchy_comb_p" %in% colnames(bt_acat$overall))
  cc <- bt_acat$overall[, "cauchy_comb_p"]
  expect_true(all(cc > 0 & cc < 1))
})

## --- 6.  The payoff: ACAT survives where d^2 is rank-degenerate -----------
##
## The reason the feature exists.  p = 40 covariates in a design with residual
## df N - B = 40 - 10 = 30 < p, so the d^2 chisquare is pinned at its rank and
## its p-value sits near 0.46 -- blind to the one strongly imbalanced covariate.
## ACAT, combining the RAW per-covariate p-values (so it must be computed before
## the Holm adjustment), still detects it.

test_that("ACAT detects imbalance where the d^2 omnibus is degenerate", {
  set.seed(11)
  K <- 10; n_s <- 4; p <- 40
  z <- blk <- numeric(0); Xl <- list()
  for (s in 1:K) {
    M <- matrix(rnorm(n_s * p), n_s, p)
    zz <- rep(0, n_s); zz[sample(n_s, 1)] <- 1
    M[zz == 1, 1] <- M[zz == 1, 1] + 4   # covariate 1 strongly imbalanced
    Xl[[s]] <- M; z <- c(z, zz); blk <- c(blk, rep(s, n_s))
  }
  X <- do.call(rbind, Xl); colnames(X) <- paste0("x", 1:p)
  dat <- data.frame(X, z = z, mset = factor(blk))
  f <- reformulate(c(paste0("x", 1:p), "strata(mset)"), response = "z")
  bt <- suppressWarnings(balanceTest(f, data = dat, cauchy.combination = TRUE))
  expect_gt(bt$overall["mset", "p.value"], 0.3)          # d^2 degenerate
  expect_lt(bt$overall["mset", "cauchy_comb_p"], 0.01)   # ACAT not fooled
})

## --- 7.  Indexing: the LAST covariate is included (bug fix) ----------------
##
## The prototype combined descriptives[-dim(descriptives)[1], ...], silently
## dropping the last row -- a real covariate when there is no missingness row.
## Put the imbalance on the last covariate; ACAT must catch it.

test_that("the last covariate is included in the combination", {
  set.seed(5); n <- 40
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    z = rep(c(1, 0), n / 2), s = factor(rep(1:4, each = n / 4)))
  dat$x3 <- ifelse(dat$z == 1, rnorm(n, 2), rnorm(n, -2))  # last covariate imbalanced
  bt <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat,
                    cauchy.combination = TRUE)
  expect_lt(bt$overall["s", "cauchy_comb_p"], 0.01)
})

## --- 8.  Missingness is a covariate: its indicator row is combined ---------
##
## RItools balances missingness as if it were another covariate, so the
## per-pattern missingness-indicator rows (named like "(x2)") carry balance
## p-values and MUST enter the combination.  Here x2 is missing far more often
## among treated units: balanced on the observed covariates, badly imbalanced on
## missingness.  ACAT over all rows catches it; dropping the missingness row (as
## the prototype's -dim()[1] indexing did) would miss it (ACAT ~ 0.6).

test_that("differential missingness is detected (missingness as covariate)", {
  set.seed(3); n <- 80
  z <- rep(c(1, 0), n / 2); s <- factor(rep(1:4, each = n / 4))
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), z = z, s = s)
  miss <- (z == 1 & runif(n) < 0.6) | (z == 0 & runif(n) < 0.05)
  d$x2[miss] <- NA
  bt <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = d,
                    cauchy.combination = TRUE)
  expect_true("(x2)" %in% rownames(bt$results))          # missingness row exists
  expect_lt(bt$overall["s", "cauchy_comb_p"], 0.01)      # and is combined in
})
