################################################################################
## Tests for the two-number magnitude balance report.
##
## Form (a): balanceMagnitude() -- a descriptive reader on top of the std.diff
##   that balanceTest() already computes.  It lives OUTSIDE the omnibus (it never
##   pre-filters covariates fed to d^2; decision D).
## Form (b): poolCRE_adjusted_percentile() -- EXPERIMENTAL, internal.  The
##   structure-adjusted whole-pool complete-randomization percentile: recover the
##   coarse confounding structure by k-means on the observed covariates,
##   residualize on it, and score the observed (matched-set adjusted) max|SMD|
##   against complete randomization on the whole pool.  Property-tested for
##   non-collapse and magnitude-sensitivity (the impossibility-pressure-test memo).
################################################################################
library("testthat")
context("magnitude balance report")

## ---- form (a): balanceMagnitude ------------------------------------------
make_bt <- function(seed = 5) {
  set.seed(seed); n <- 60
  d <- data.frame(z = rep(c(1, 0), 30), s = factor(rep(1:30, each = 2)),
                  x1 = rnorm(n), x2 = rnorm(n) + rep(c(0.4, 0), 30), x3 = rnorm(n))
  balanceTest(z ~ x1 + x2 + x3 + strata(s), data = d)
}

test_that("balanceMagnitude exists and reads std.diff from $results", {
  expect_true(is.function(balanceMagnitude))
  bt <- make_bt()
  bm <- balanceMagnitude(bt, threshold = 0.25)
  expect_s3_class(bm, "balanceMagnitude")
  sd_s <- bt$results[, "std.diff", "s"]
  tab  <- bm[["s"]]
  expect_setequal(round(tab$std.diff, 8), round(as.numeric(sd_s), 8))
  expect_equal(tab$abs.std.diff, abs(tab$std.diff))
})

test_that("balanceMagnitude reports the global max|SMD| and flips exceedances on the threshold", {
  bt <- make_bt()
  sd_s <- bt$results[, "std.diff", "s"]
  mn <- min(abs(sd_s)); mx <- max(abs(sd_s))
  tab  <- balanceMagnitude(bt)[["s"]]
  expect_equal(attr(tab, "max.abs.std.diff"), mx)
  ## below the smallest |SMD| -> all exceed; above the largest -> none exceed
  expect_true(all(balanceMagnitude(bt, threshold = mn * 0.5)[["s"]]$exceeds))
  expect_false(any(balanceMagnitude(bt, threshold = mx * 1.5)[["s"]]$exceeds))
})

test_that("balanceMagnitude covers every stratification including the unstratified --", {
  bt <- make_bt()
  bm <- balanceMagnitude(bt)
  expect_true(all(c("s", "--") %in% names(bm)))
  ## a specific stratification can be requested
  expect_identical(names(balanceMagnitude(bt, stratification = "--")), "--")
  expect_error(balanceMagnitude(bt, stratification = "nope"), "unknown")
})

test_that("balanceMagnitude is in the package exports", {
  expect_true("balanceMagnitude" %in% getNamespaceExports("RItools"))
})

## ---- form (b): poolCRE_adjusted_percentile (EXPERIMENTAL, internal) -------
## self-contained DGP (tests cannot source vignettes/): G coarse groups x sets x
## (2T+2C); fixed per-covariate gap g along u; within-set noise scaled by shrink.
mk_matched <- function(shrink, g = 0.2, K = 3, group_centers = c(0, 8, 16, 24),
                       per = 5, seed = 1) {
  set.seed(seed)
  u <- rep(1, K) / sqrt(K); Xl <- list(); ms <- grp <- z <- integer(0); sid <- 0L
  for (gi in seq_along(group_centers)) for (s in seq_len(per)) {
    sid <- sid + 1L; ctr <- group_centers[gi] + rnorm(K) * 1.5
    blk <- rbind(ctr + (g / 2) * u + shrink * rnorm(K),
                 ctr + (g / 2) * u + shrink * rnorm(K),
                 ctr - (g / 2) * u + shrink * rnorm(K),
                 ctr - (g / 2) * u + shrink * rnorm(K))
    Xl[[sid]] <- blk; ms <- c(ms, rep(sid, 4)); grp <- c(grp, rep(gi, 4)); z <- c(z, c(1, 1, 0, 0))
  }
  X <- do.call(rbind, Xl); colnames(X) <- paste0("x", seq_len(K))
  list(X = X, z = z, set = factor(ms), group = factor(grp))
}

test_that("poolCRE_adjusted_percentile does NOT collapse as the match tightens", {
  skip_on_cran()
  f <- RItools:::poolCRE_adjusted_percentile
  dloose <- mk_matched(0.60, seed = 11); dtight <- mk_matched(0.03, seed = 11)
  set.seed(1); p_loose <- f(dloose$X, dloose$z, dloose$set, k = 4, B = 400)
  set.seed(1); p_tight <- f(dtight$X, dtight$z, dtight$set, k = 4, B = 400)
  ## the design beats a whole-pool coin flip at both tightness levels; crucially
  ## the percentile does NOT climb toward 1 as the match tightens (the omnibus p would)
  expect_lt(p_tight, 0.9)
  expect_lt(p_tight - p_loose, 0.5)
})

test_that("poolCRE_adjusted_percentile is magnitude-sensitive (rises with a real gap)", {
  skip_on_cran()
  f <- RItools:::poolCRE_adjusted_percentile
  d0 <- mk_matched(0.10, g = 0.0, seed = 7); d1 <- mk_matched(0.10, g = 1.0, seed = 7)
  set.seed(2); p0 <- f(d0$X, d0$z, d0$set, k = 4, B = 400)
  set.seed(2); p1 <- f(d1$X, d1$z, d1$set, k = 4, B = 400)
  expect_gt(p1, p0 + 0.1)
})

test_that("k-means structure recovery roughly matches the oracle group adjustment", {
  skip_on_cran()
  f <- RItools:::poolCRE_adjusted_percentile
  d <- mk_matched(0.10, g = 0.6, seed = 3)
  ## oracle: residualize on the TRUE group label, same statistic/reference
  Xr_oracle <- apply(d$X, 2, function(x) x - ave(x, d$group))
  set.seed(4); p_km <- f(d$X, d$z, d$set, k = 4, B = 400)
  ## recompute the oracle percentile inline with the same internals
  sdp <- RItools:::.magnitude_pooled_sd(Xr_oracle, d$z)
  Tobs <- RItools:::.within_stratum_max_smd(Xr_oracle, d$set, d$z, sdp)
  pool <- factor(rep(1L, length(d$z)))
  set.seed(4); Tnull <- replicate(400, RItools:::.within_stratum_max_smd(Xr_oracle, pool, sample(d$z), sdp))
  p_oracle <- mean(Tnull <= Tobs - 1e-12)
  expect_lt(abs(p_km - p_oracle), 0.3)
})

## ---- form (b): data-driven k and selectable centering (Jake's decisions) --

test_that(".choose_k_gap recovers ~the true number of coarse groups", {
  skip_on_cran()
  ## 4 well-separated coarse groups -> the gap statistic should pick about 4
  d <- mk_matched(0.10, g = 0.6, seed = 3)
  set.seed(9)
  ksel <- RItools:::.choose_k_gap(d$X)
  expect_true(ksel >= 3 && ksel <= 5)
})

test_that("data-driven k (k = NULL, the default) stays non-collapsing and magnitude-sensitive", {
  skip_on_cran()
  f <- RItools:::poolCRE_adjusted_percentile
  d0 <- mk_matched(0.10, g = 0.0, seed = 7); d1 <- mk_matched(0.10, g = 1.0, seed = 7)
  set.seed(2); p0 <- f(d0$X, d0$z, d0$set, k = NULL, B = 300)
  set.seed(2); p1 <- f(d1$X, d1$z, d1$set, k = NULL, B = 300)
  expect_gt(p1, p0 + 0.1)
})

test_that("centering = 'uniform' is a valid (uniform) percentile under complete randomization", {
  skip_on_cran()
  ## when the observed assignment is itself a complete-randomization draw and the
  ## reference is complete randomization, the pool/pool percentile is ~Uniform
  f <- RItools:::poolCRE_adjusted_percentile
  d <- mk_matched(0.30, g = 0.0, seed = 20)            # null: no systematic gap
  set.seed(5)
  ps <- replicate(30, f(d$X, sample(d$z), d$set, k = 4, centering = "uniform", B = 200))
  expect_equal(mean(ps), 0.5, tolerance = 0.12)
})

test_that("centering = 'conservative' is downward-biased relative to 'uniform'", {
  skip_on_cran()
  ## the matched-set-centered observed statistic is smaller than the pool-centered
  ## reference, so under a complete-randomization observed assignment the
  ## conservative percentile sits below the uniform one (paired on the same draw)
  f <- RItools:::poolCRE_adjusted_percentile
  d <- mk_matched(0.30, g = 0.0, seed = 22)
  set.seed(7)
  diffs <- replicate(20, {
    zz <- sample(d$z)                                  # complete-randomization draw
    pu <- f(d$X, zz, d$set, k = 4, centering = "uniform",      B = 200)
    pc <- f(d$X, zz, d$set, k = 4, centering = "conservative", B = 200)
    pu - pc
  })
  expect_gt(mean(diffs), 0)
})
