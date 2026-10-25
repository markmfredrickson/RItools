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
## That conversion sends p_k = 1 to minus infinity, so one covariate whose
## adjusted difference is exactly 0 (p_k = 1) sets p_ACAT to 1 whatever the other
## covariates say.  acat_pvalue() therefore uses, by default, the truncated
## conversion of Gui, Jiang, and Wang (2025, Biometrika 112(4), asaf038), as
## riposte does: with truncation share t (0.9 by default) each p_k becomes
## tan((0.5 - t p_k) pi), so p_k = 1 becomes tan(-0.4 pi) = -3.08, and with n
## converted values summing to S the combined p-value is min(1, n P(C > S) / t)
## for a standard Cauchy C.  truncation = 1 gives Liu and Xie's version.  The
## cost, measured in vignettes/truncated-cauchy-memo.md: at the 0.05 level the
## truncated version rejects a true null somewhat more often (0.057 with 6
## independent covariates); at 0.01 it stays near 0.01.
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

test_that("Liu and Xie's combination reproduces an identical p-value (idempotence)", {
  acat_lx <- function(p) RItools:::acat_pvalue(p, truncation = 1)
  expect_equal(acat_lx(rep(0.30, 5)), 0.30, tolerance = 1e-8)
  expect_equal(acat_lx(rep(0.73, 3)), 0.73, tolerance = 1e-8)
  expect_equal(acat_lx(rep(0.05, 8)), 0.05, tolerance = 1e-8)
})

## The truncated combination is idempotent only for a single p-value: with
## n = 1, P(C > tan((0.5 - t p) pi)) = t p, so the combined p-value is
## t p / t = p for any truncation share t.  (Five copies of 0.3 combine to
## 0.394, not 0.3; see vignettes/truncated-cauchy-memo.md.)

test_that("the truncated combination returns a single p-value unchanged", {
  for (p in c(1e-200, 0.003, 0.2, 0.7, 1)) {
    expect_equal(RItools:::acat_pvalue(p), p, tolerance = 1e-12)
    expect_equal(RItools:::acat_pvalue(p, truncation = 0.5), p, tolerance = 1e-12)
  }
})

## --- 2.  Null calibration: a combined p-value must be uniform --------------
##
## Under the global null every per-covariate p_k is Uniform(0,1).  A valid
## combination is then itself Uniform(0,1): mean 0.5 and rejection rates equal to
## their nominal levels.  This is the property that justifies using the number as
## a p-value at all, so it is the most important test here.

test_that("Liu and Xie's combination is calibrated under the global null", {
  skip_on_cran()                       # Monte Carlo; deterministic via the seed
  set.seed(2024)
  B <- 4000; K <- 6
  null_p <- replicate(B, RItools:::acat_pvalue(runif(K), truncation = 1))
  expect_equal(mean(null_p), 0.5, tolerance = 0.02)
  expect_equal(mean(null_p < 0.05), 0.05, tolerance = 0.015)
  expect_equal(mean(null_p < 0.10), 0.10, tolerance = 0.02)
  expect_gt(suppressWarnings(ks.test(null_p, "punif")$p.value), 0.01)
})

## The truncated combination is not uniform under the null: min(1, ...) puts
## mass at 1, and the heavy-tail approximation behind n P(C > S) / t is
## accurate only for small p-values.  What it must deliver is the level at small
## alpha.  With 20,000 draws the standard error of a rejection rate near 0.01 is
## sqrt(0.01 * 0.99 / 20000) = 0.0007, so we accept 0.01 +/- 0.003.  At 0.05 the
## memo measured 0.057 with K = 6 independent covariates; the bound 0.0625
## records that accepted inflation and would catch a larger one.

test_that("the truncated combination keeps its level at 0.01 under the global null", {
  skip_on_cran()                       # Monte Carlo; deterministic via the seed
  set.seed(2025)
  B <- 20000; K <- 6
  null_p <- replicate(B, RItools:::acat_pvalue(runif(K)))
  expect_equal(mean(null_p < 0.01), 0.01, tolerance = 0.003)
  expect_lt(mean(null_p < 0.05), 0.0625)
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

test_that("acat_pvalue handles boundaries, NAs, and stays in (0,1]", {
  expect_true(is.finite(RItools:::acat_pvalue(c(0, 0.5, 1))))
  ## NAs are dropped: the one remaining p-value combines to itself
  expect_equal(RItools:::acat_pvalue(c(0.3, NA)), 0.30, tolerance = 1e-8)
  expect_equal(RItools:::acat_pvalue(c(0.3, NA, 0.3), truncation = 1), 0.30,
               tolerance = 1e-8)
  expect_true(is.na(RItools:::acat_pvalue(c(NA, NA))))
  set.seed(7)
  ## the truncated combination can return exactly 1 (from its min(1, ...))
  for (i in 1:50) {
    val <- RItools:::acat_pvalue(runif(sample(2:8, 1)))
    expect_true(val > 0 && val <= 1)
  }
  expect_error(RItools:::acat_pvalue(0.5, truncation = 0), "truncation")
  expect_error(RItools:::acat_pvalue(0.5, truncation = 1.1), "truncation")
})

## --- 4b.  The truncation: a p-value of 1 no longer decides the result -----
##
## Four covariates with p = 0.01 and one with p = 1.  Liu and Xie's conversion
## sends the 1 to minus infinity (capped at about -1.4e15), so the combination
## is about 1.  The truncated conversion sends it to tan(-0.4 pi) = -3.08, and
## the combined p-value follows the formula min(1, n P(C > S) / t).

test_that("one p-value of 1 does not set the truncated combination to 1", {
  p <- c(0.01, 0.01, 0.01, 0.01, 1)
  expect_gt(RItools:::acat_pvalue(p, truncation = 1), 0.999)
  S <- 4 * tan((0.5 - 0.9 * 0.01) * pi) + tan(-0.4 * pi)
  expect_equal(RItools:::acat_pvalue(p), 5 * pcauchy(S, lower.tail = FALSE) / 0.9)
  expect_lt(RItools:::acat_pvalue(p), 0.05)
})

## A p-value too small for tan() to convert (1e-300 converts to about
## 1/(0.9e-300 pi), which overflows) must still give a representable result:
## one tiny p-value among n gives about n times that p-value.

test_that("a tiny p-value neither overflows nor is lost", {
  res <- RItools:::acat_pvalue(c(1e-300, 0.5, 0.5))
  expect_gt(res, 0)
  expect_equal(res / 3e-300, 1, tolerance = 1e-6)
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
  expect_true(all(cc > 0 & cc <= 1))   # the truncated version can return exactly 1
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

## --- 9.  ACAT combines raw p-values, whatever p.adjust.method is ----------
##
## The per-covariate p-values reported in $results are multiplicity-adjusted
## (Holm by default).  ACAT must combine the unadjusted ones: Holm sets many
## p-values to exactly 1, and each 1 drives the Cauchy sum toward minus
## infinity.  If the combination read the adjusted column, changing
## p.adjust.method would change cauchy_comb_p.  This guards the ordering of
## the ACAT block relative to the p.adjust() loop in balanceTest().

test_that("cauchy_comb_p does not depend on p.adjust.method", {
  set.seed(5); n <- 40
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    z = rep(c(1, 0), n / 2), s = factor(rep(1:4, each = n / 4)))
  dat$x3 <- ifelse(dat$z == 1, rnorm(n, 1), rnorm(n, -1))
  bt_holm <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat,
                         cauchy.combination = TRUE, p.adjust.method = "holm")
  bt_none <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat,
                         cauchy.combination = TRUE, p.adjust.method = "none")
  expect_equal(bt_holm$overall[, "cauchy_comb_p"],
               bt_none$overall[, "cauchy_comb_p"])
  ## and the combination equals ACAT applied to the unadjusted p-values
  expect_equal(unname(bt_none$overall["s", "cauchy_comb_p"]),
               RItools:::acat_pvalue(bt_none$results[, "p", "s"]))
})

## --- 10.  Integration: an exactly balanced covariate does not hide imbalance
##
## xb is binary with the same treated and control share (3/5) in every block,
## so its adjusted difference is exactly 0 and its p-value is exactly 1.  x1 is
## imbalanced (p = 1.5e-05).  Under Liu and Xie's conversion cauchy_comb_p was
## 0.99999999999999933; the truncated conversion gives 4.6e-05.

test_that("an exactly balanced covariate does not hide an imbalanced one", {
  set.seed(2); n <- 40
  d2 <- data.frame(z = rep(c(1, 0), n / 2), s = factor(rep(1:4, each = n / 4)))
  d2$xb <- rep(c(1, 1, 0, 0), n / 4)
  d2$x1 <- ifelse(d2$z == 1, rnorm(n, 1), rnorm(n, -1))
  d2$x2 <- rnorm(n)
  bt <- balanceTest(z ~ xb + x1 + x2 + strata(s), data = d2,
                    cauchy.combination = TRUE)
  expect_equal(unname(bt$results["xb", "p", "s"]), 1)
  expect_lt(bt$overall["s", "cauchy_comb_p"], 1e-3)
})
