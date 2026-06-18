################################################################################
## Tests-first for the ACAT-omnibus switch guidance.
##
## ACAT itself -- the cauchy.combination = TRUE path: the published formula, null
## calibration, single-signal sensitivity, high-dimensional detection, and
## opt-in/off-by-default behavior -- is tested in test.cauchycomb.R.  THIS file
## tests only the NEW behavior we are adding:
##
##   when the d^2 omnibus is degenerate or near-degenerate (#covariates >=
##   within-block residual df = rank(V_d)), where the chi-square statistic
##   freezes at its rank and its p-value parks near 0.45, balanceTest() should
##   WARN the user to switch to the ACAT omnibus (cauchy.combination = TRUE).
##
## Today the only warning in that regime is the generic "Degrees of freedom
## exceeds units less number of strata ..." (R/Design.R check_for_degenerate),
## which never mentions ACAT.  So the ACAT-suggestion (test 1) is expected RED;
## that SOME degeneracy warning fires (test 2) and that ACAT stays informative
## there (test 3) are GREEN and document why the suggestion is worth making.
################################################################################
library("testthat")
context("ACAT omnibus switch guidance")

## High-dimensional matched sets: B blocks of size block_size, K covariates.  The
## within-block residual df is N - B = B*(block_size - 1); with K above that the
## d^2 covariance is rank-deficient.  A real confound sits on covariate
## `confound_col` (treated shifted up by `confound`).
make_highdim <- function(B, block_size, K, confound_col = 1, confound = 4, seed = 11) {
  set.seed(seed)
  z <- blk <- numeric(0); Xl <- list()
  for (s in seq_len(B)) {
    M  <- matrix(rnorm(block_size * K), block_size, K)
    zz <- rep(0, block_size); zz[sample(block_size, 1)] <- 1
    M[zz == 1, confound_col] <- M[zz == 1, confound_col] + confound   # real imbalance
    Xl[[s]] <- M; z <- c(z, zz); blk <- c(blk, rep(s, block_size))
  }
  X <- do.call(rbind, Xl); colnames(X) <- paste0("x", seq_len(K))
  data.frame(X, z = z, mset = factor(blk))
}
hd_form <- function(K) reformulate(c(paste0("x", 1:K), "strata(mset)"), response = "z")

## --- 1. the degenerate-omnibus warning should point users to ACAT -----------
## 20 covariates, 15 pairs: residual df = 15 < 20, so d^2 is rank-degenerate.
## The feature: the warning names cauchy.combination / ACAT so the user knows the
## informative alternative.  RED today (the only warning is the generic df one).
test_that("a degenerate omnibus warns the user to consider the ACAT omnibus", {
  d <- make_highdim(B = 15, block_size = 2, K = 20)
  expect_warning(balanceTest(hd_form(20), data = d),
                 regexp = "cauchy|ACAT|Cauchy")
})

## --- 2. (guard, GREEN) a degenerate omnibus warns at all --------------------
## Documents the current behavior the new message builds on: SOME warning fires.
test_that("a degenerate omnibus warns about degrees of freedom (current behavior)", {
  d <- make_highdim(B = 15, block_size = 2, K = 20)
  expect_warning(balanceTest(hd_form(20), data = d),
                 regexp = "[Dd]egrees of freedom")
})

## --- 3. (GREEN) the reason for the suggestion: ACAT stays informative --------
## With a strong confound on x1, the frozen d^2 cannot reject (p ~ 0.45) but the
## ACAT combination of the per-covariate p-values does.  This is the payoff that
## justifies steering users to ACAT when the d^2 omnibus is degenerate.
test_that("ACAT rejects a real imbalance where the degenerate d^2 cannot", {
  d  <- make_highdim(B = 15, block_size = 2, K = 20)
  bt <- suppressWarnings(balanceTest(hd_form(20), data = d, cauchy.combination = TRUE))
  expect_gt(bt$overall["mset", "p.value"], 0.3)          # d^2 frozen near its rank
  expect_lt(bt$overall["mset", "cauchy_comb_p"], 0.05)   # ACAT not fooled
})

## --- 4. a well-conditioned omnibus does NOT warn about ACAT -----------------
## The suggestion must be specific to the degenerate regime: a low-dimensional,
## full-rank problem should not nag users to switch.  (Today there is no ACAT
## warning at all, so this passes vacuously now; it guards against a future
## over-eager warning.)
test_that("a well-conditioned omnibus does not suggest ACAT", {
  d <- make_highdim(B = 30, block_size = 2, K = 3)       # residual df 30 >> K = 3
  expect_no_warning(balanceTest(hd_form(3), data = d))
})
