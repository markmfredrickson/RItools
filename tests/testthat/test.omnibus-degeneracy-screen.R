################################################################################
## Tests-first for the numerical-stability screen on the d^2 omnibus.
##
## Substantive principle.  A covariate (or direction) with negligible
## WITHIN-STRATUM variance is, by construction of the matching, already balanced:
## there is nothing for the within-stratum permutation test to test on it, and
## inverting its ~0 null variance is what makes the omnibus numerically
## unstable / undefined.  balanceTest() should therefore
##   (1) DROP such a covariate from the omnibus, lowering the reported df, and
##       MESSAGE the user, naming the covariate, e.g.
##       "within strata there is too little variance in x1 for balance testing
##        -- it is nearly perfectly balanced";
##   (2) DEGRADE GRACEFULLY (no error) when EVERY covariate is negligible,
##       returning a result that signals "no testable covariates";
##   (3) make the screening decision from the DESIGN only (covariate values +
##       strata), so it is invariant to the treatment assignment (ungameable);
##   (4) leave a clean (non-degenerate) problem's omnibus UNCHANGED.
##
## Current code (R/Design.R ~991) screens only ssvar <= .Machine$double.eps -- an
## ABSOLUTE machine-epsilon cutoff, with no message -- and R/utils.R
## XtX_pseudoinv_sqrt() ERRORS (stop(), ~line 381) when nothing survives.  So:
##   - the MESSAGE (1) and the RELATIVE tolerance (1b) are expected RED;
##   - GRACEFUL DEGRADATION (2) is expected RED (currently errors);
##   - assignment-invariance (3) and non-interference (4) already hold for
##     matched pairs and are GREEN regression guards for the new screen.
################################################################################
library("testthat")
context("omnibus numerical-stability screen")

## Matched-pairs builder.  S pairs, K covariates.  within_sd[j] is the
## within-pair SD of covariate j (small => nearly exactly matched); pair
## locations are N(0,1) -- the between-pair variation that stratification removes,
## so each covariate's TOTAL variance stays ~1 even when its within-pair variance
## is tiny.  `confound` adds a constant treated-minus-control shift to x1.
make_pairs_df <- function(S, within_sd, confound = 0, seed = 1) {
  set.seed(seed)
  K   <- length(within_sd)
  loc <- matrix(rnorm(S * K), S, K)                          # pair locations
  dif <- matrix(rnorm(S * K), S, K) *
         matrix(within_sd, S, K, byrow = TRUE)               # within-pair differences
  dif[, 1] <- dif[, 1] + confound                            # systematic tilt on x1
  X <- matrix(0, 2 * S, K)
  X[seq(1, 2 * S, 2), ] <- loc + dif / 2                     # treated rows
  X[seq(2, 2 * S, 2), ] <- loc - dif / 2                     # control rows
  out <- data.frame(z = rep(c(1L, 0L), S), s = factor(rep(seq_len(S), each = 2)))
  for (j in seq_len(K)) out[[paste0("x", j)]] <- X[, j]
  out
}
bt_form <- function(K) reformulate(c(paste0("x", 1:K), "strata(s)"), response = "z")
## covariates excluded from the omnibus show up as z = NA in $results
screened_vars <- function(bt, strat = "s") names(which(is.na(bt$results[, "z", strat])))

## --- 1a. an EXACTLY-matched covariate is dropped, WITH a message -------------
## x1 identical within every pair (within_sd 0); x2, x3 vary.  The current code
## already EXCLUDES x1 (its null variance is exactly 0), so df = 2 and x1's z is
## NA are GREEN -- but it does so SILENTLY.  The message naming x1 is the RED.
test_that("an exactly-matched covariate is screened out, and the user is told", {
  d <- make_pairs_df(S = 40, within_sd = c(0, 1, 1))
  expect_message(bt <- balanceTest(bt_form(3), data = d),
                 regexp = "x1")                       # RED: no message today
  expect_equal(unname(bt$overall["s", "df"]), 2)      # GREEN: x1 already excluded
  expect_true("x1" %in% screened_vars(bt))            # GREEN
})

## --- 1b. a RELATIVELY negligible covariate is dropped (relative tolerance) ---
## x1's within-pair SD is 1e-5: its within-stratum variance (~1e-10) is negligible
## next to its unit total variance, but well ABOVE machine epsilon -- so the
## current absolute cutoff does NOT catch it and x1 stays in the omnibus (df = 3).
## A relative tolerance should drop it (df 3 -> 2) and message.  Both are RED.
test_that("a relatively-negligible covariate is screened out (relative tolerance)", {
  d <- make_pairs_df(S = 40, within_sd = c(1e-5, 1, 1))
  expect_message(bt <- balanceTest(bt_form(3), data = d),
                 regexp = "x1")                       # RED
  expect_equal(unname(bt$overall["s", "df"]), 2)      # RED: today df = 3
  expect_true("x1" %in% screened_vars(bt))            # RED
})

## --- 2. graceful degradation: NO error when every covariate is degenerate ----
## All three covariates have ~0 within-pair variance (within_sd 1e-10), so the
## surviving covariance is empty and XtX_pseudoinv_sqrt() currently stop()s.
## balanceTest() must not error; it should signal "no testable covariates".  We
## catch the error so this reports as a clean FAILURE rather than an ERROR.
test_that("balanceTest does not error when all covariates are (near) exactly matched", {
  d <- make_pairs_df(S = 30, within_sd = rep(1e-10, 3))
  res <- tryCatch(suppressMessages(balanceTest(bt_form(3), data = d)),
                  error = function(e) structure(list(msg = conditionMessage(e)),
                                                class = "caught_error"))
  expect_false(inherits(res, "caught_error"))         # RED: today it stop()s
})

## --- 2b. graceful abstention when all covariates are relatively negligible ---
## within_sd 1e-5 on all three: above machine eps, so no error today (the omnibus
## is computed on pure noise, df = 3).  The feature should instead abstain --
## report no testable covariates (df 0 or NA p.value) -- and say so.
test_that("balanceTest abstains (does not test noise) when all covariates are negligible", {
  d <- make_pairs_df(S = 30, within_sd = rep(1e-5, 3))
  expect_message(bt <- balanceTest(bt_form(3), data = d),
                 regexp = "balance")                  # RED: no message today
  abstained <- is.na(bt$overall["s", "p.value"]) || unname(bt$overall["s", "df"]) == 0
  expect_true(abstained)                              # RED: today df = 3, p finite
})

## --- 3. (guard, GREEN) the screen is assignment-invariant --------------------
## The set of screened covariates depends on covariate values and strata, NOT on
## the treatment vector.  For matched pairs this already holds (the per-stratum
## n1 = n0 = 1 regardless of which unit is treated).  Flip treatment within every
## pair and the screened set must be unchanged -- a property the new relative
## screen must also keep, which is what makes it ungameable.
test_that("the screened-covariate set is invariant to within-pair re-randomization", {
  d <- make_pairs_df(S = 30, within_sd = c(1e-10, 1, 1), confound = 0.4)
  scr_obs <- suppressMessages(screened_vars(balanceTest(bt_form(3), data = d)))
  d_flip <- d
  flip <- rep(rbinom(nlevels(d$s), 1, 0.5) == 1, each = 2)
  d_flip$z <- ifelse(flip, 1L - d_flip$z, d_flip$z)
  scr_flip <- suppressMessages(screened_vars(balanceTest(bt_form(3), data = d_flip)))
  expect_identical(sort(scr_obs), sort(scr_flip))
})

## --- 4. (guard, GREEN) non-interference on a clean problem -------------------
## With no degenerate covariate the omnibus must be the full-rank chi-square on
## all K covariates -- the screen must not touch it.
test_that("a non-degenerate omnibus is left intact (df = number of covariates)", {
  d <- make_pairs_df(S = 50, within_sd = c(1, 1, 1), confound = 0.3)
  bt <- balanceTest(bt_form(3), data = d)
  expect_equal(unname(bt$overall["s", "df"]), 3)
  expect_true(is.finite(bt$overall["s", "chisquare"]) && bt$overall["s", "chisquare"] > 0)
  expect_true(bt$overall["s", "p.value"] > 0 && bt$overall["s", "p.value"] < 1)
  expect_length(screened_vars(bt), 0)
})
