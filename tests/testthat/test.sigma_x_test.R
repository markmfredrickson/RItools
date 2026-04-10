################################################################################
## Tests for the sigma_x balance test:  T = d' Sigma_x^{-1} d
##
## The new omnibus statistic standardizes the adjusted differences vector by an
## external (or sample) covariance Sigma_x rather than by the permutation
## covariance Cov(d).  Under the strata-conditional randomization null,
## T ~ sum_k lambda_k * chi^2_1 with lambda_k = eigenvalues of
## Sigma_x^{-1/2} Cov(d) Sigma_x^{-1/2}, but only ASYMPTOTICALLY (when d is
## approximately MVN).  At finite n, the second-order moment of T can deviate
## substantially from 2 * sum_k lambda_k^2; see Section 6.
##
## Tests use exact within-stratum enumeration on tiny examples as ground truth.
## Helpers are in helper-sigma_x_test.R.
################################################################################
library("testthat")
context("sigma_x balance test")

## --- 1.  pvalue_quad_form / sigma_x_pvalue --------------------------------
##
## Lowest-level: given the test stat value q and either lambda or (E[T],
## Var[T]), return P(T > q) under the chosen null backend.

test_that("sigma_x_pvalue 'imhof' backend matches CompQuadForm::imhof", {
  ## Anchor: our wrapper must delegate correctly to CompQuadForm.  We are not
  ## re-implementing imhof; we just want to be sure we hand the eigenvalues in
  ## with the right sign convention and pull the right field out of the result.
  skip_if_not_installed("CompQuadForm")
  lambda <- c(2, 1, 0.5)
  q <- 4
  expect_equal(
    sigma_x_pvalue(q, lambda = lambda, method = "imhof"),
    CompQuadForm::imhof(q, lambda = lambda)$Qq,
    tolerance = 1e-6
  )
})

test_that("sigma_x_pvalue 'satterthwaite_asymptotic' uses (sum lambda, 2 sum lambda^2)", {
  ## Asymptotic Satterthwaite approximates T ~ a chi^2_v with
  ##   E[T] = sum lambda,  Var[T] = 2 sum lambda^2,
  ##   v = (sum lambda)^2 / (sum lambda^2),  a = (sum lambda^2) / (sum lambda).
  ## When all eigenvalues are equal to c, a = c and v = p, so T ~ c chi^2_p
  ## exactly under Gaussian d -- this is the cleanest sanity check.
  lambda <- rep(2, 5)  # five equal eigenvalues
  q <- 12
  expect_equal(
    sigma_x_pvalue(q, lambda = lambda, method = "satterthwaite_asymptotic"),
    pchisq(q / 2, df = 5, lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("sigma_x_pvalue 'satterthwaite_finite' uses supplied (E[T], Var[T])", {
  ## Finite-sample backend: caller supplies E[T] and Var[T] computed under the
  ## EXACT randomization distribution (not via the Gaussian identity).  We
  ## moment-match a chi-square: a = Var/(2 E),  v = 2 E^2 / Var,  T ~ a chi^2_v.
  ## Hand check with a known case.
  EM   <- 6
  VarM <- 12        # so a = 1, v = 6, T ~ chi^2_6 exactly
  q    <- 4
  expect_equal(
    sigma_x_pvalue(q, EM = EM, VarM = VarM, method = "satterthwaite_finite"),
    pchisq(q, df = 6, lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("sigma_x_pvalue 'imhof' and 'satterthwaite_asymptotic' agree on equal eigenvalues", {
  ## When eigenvalues are equal, Imhof and the asymptotic Satterthwaite must
  ## give identical p-values: T ~ c chi^2_p in both cases.
  skip_if_not_installed("CompQuadForm")
  lambda <- rep(1.5, 4)
  for (q in c(0.1, 1, 5, 12, 30)) {
    expect_equal(
      sigma_x_pvalue(q, lambda = lambda, method = "imhof"),
      sigma_x_pvalue(q, lambda = lambda, method = "satterthwaite_asymptotic"),
      tolerance = 1e-5,
      info = paste("q =", q)
    )
  }
})

## --- 2.  Closed form for Cov(d_raw) ---------------------------------------
##
## Section 3.1 of the design.  Cov(d_raw) computed as
##   sum_s [ n1_s n0_s / (n_s (n_s - 1)) ] * S_xs
## must equal the empirical covariance over exact enumeration to floating-
## point precision.  This is verified inside the test rather than relying on
## the helper to be correct -- the helper is the same formula, so we test the
## *package's* implementation, not the helper's.

test_that("randomization_cov_d closed form matches exact enumeration (equal strata)", {
  fx <- make_small_fixture()
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  emp <- exact_randomization_dist(zs, fx$X, fx$strata,
                                  sigma_x = diag(fx$p))$V_d_emp
  pkg <- randomization_cov_d(fx$X, fx$strata, fx$n1_per_stratum)
  expect_equal(pkg, emp, tolerance = 1e-12, check.attributes = FALSE)
})

test_that("randomization_cov_d closed form matches exact enumeration (uneven strata)", {
  fx <- make_uneven_fixture()
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  emp <- exact_randomization_dist(zs, fx$X, fx$strata,
                                  sigma_x = diag(fx$p))$V_d_emp
  pkg <- randomization_cov_d(fx$X, fx$strata, fx$n1_per_stratum)
  expect_equal(pkg, emp, tolerance = 1e-12, check.attributes = FALSE)
})

test_that("randomization_cov_d agrees with HB08 tcov (up to scaling) in main", {
  ## HB08's tcov computes Cov(d) for the unnormalized ssn = (z - pi)' x_tilde
  ## vector in x_tilde units.  Our randomization_cov_d computes Cov(d_raw) in
  ## raw covariate units.  The two are related by the StrataWeightRatio.  We
  ## test in a *trivial* equal-weight case where the ratio is constant, so the
  ## two coincide up to a known scalar.
  data(nuclearplants, package = "RItools")
  bt <- balanceTest(pr ~ date + t1 + t2, data = nuclearplants)
  hb_tcov <- attr(bt$overall, "tcov")[["--"]]

  ## Unstratified: closed form V_d for raw d = (z - pi)' X
  X <- as.matrix(nuclearplants[, c("date", "t1", "t2")])
  strata <- factor(rep("--", nrow(nuclearplants)))
  n1 <- sum(nuclearplants$pr)
  V_d <- randomization_cov_d(X, strata, n1_per_stratum = n1)

  ## hb_tcov is in HB08 internal units; it also includes one all-zero
  ## row/column for the NotMissing indicator (which is all-ones, hence
  ## centered to 0) when there is no actual missingness, so it can have a
  ## strictly larger column count than pkg_V_d.  The eigenvalues of the two
  ## matrices, once each is normalized by its trace and zero eigenvalues are
  ## dropped, must match (both are valid covariance matrices for the same
  ## Gaussian limit, differing only by an overall positive scalar that cancels
  ## in the eigen-direction sense).
  drop_zero <- function(x, tol = sqrt(.Machine$double.eps)) {
    sort(x[x > tol * max(x, 1)])
  }
  hb_eig  <- drop_zero(eigen(hb_tcov / sum(diag(hb_tcov)), symmetric = TRUE,
                             only.values = TRUE)$values)
  pkg_eig <- drop_zero(eigen(V_d / sum(diag(V_d)), symmetric = TRUE,
                             only.values = TRUE)$values)
  expect_equal(pkg_eig, hb_eig, tolerance = 1e-8)
})

## --- 3.  Default Sigma_x and user supplied --------------------------------

test_that("default_sigma_x equals (1/(N-K)) * sum_s SSE within strata", {
  fx <- make_small_fixture()
  expected <- within_stratum_pooled_cov_helper(fx$X, fx$strata)
  got      <- default_sigma_x(fx$X, fx$strata)
  expect_equal(got, expected, tolerance = 1e-12, check.attributes = FALSE)
})

test_that("default_sigma_x reduces to cov() in the unstratified case", {
  ## With one stratum, within-stratum-pooled with N - K = N - 1 in the
  ## denominator must equal the ordinary unbiased sample covariance.
  set.seed(7); n <- 30; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  strata <- factor(rep("a", n))
  expect_equal(default_sigma_x(X, strata), cov(X),
               tolerance = 1e-12, check.attributes = FALSE)
})

test_that("user-supplied sigma_x is honored byte-for-byte", {
  fx <- make_small_fixture()
  ## Pick an obviously not-the-default Sigma_x: 5 * I.
  custom <- 5 * diag(fx$p)
  z <- numeric(fx$n); z[c(1, 3, 5, 7)] <- 1  # one assignment

  res <- sigma_x_test(fx$X, fx$strata, z,
                      sigma_x = custom,
                      null = "satterthwaite_asymptotic")
  ## Read the sigma_x out of the result; must match the supplied matrix
  expect_equal(res$sigma_x_used, custom, tolerance = 1e-12,
               check.attributes = FALSE)

  ## And the test stat itself must equal d' (5 I)^{-1} d = ||d||^2 / 5
  pi_i <- ave(z, fx$strata)
  d <- drop((z - pi_i) %*% fx$X)
  expect_equal(res$statistic, sum(d^2) / 5, tolerance = 1e-12)
})

## --- 4.  Test stat math: invariance, reductions ---------------------------

test_that("sigma_x_test reduces to HB08 d^2 when sigma_x = Cov(d)", {
  ## Mathematically: with Sigma_x = V_d, Sigma_x^{-1/2} V_d Sigma_x^{-1/2} = I,
  ## all eigenvalues equal 1, and T = d' V_d^{-1} d, which is exactly the
  ## Hansen-Bowers d^2 statistic.  This must equal HB08's chisquare on the same
  ## data.
  fx <- make_small_fixture()
  z  <- numeric(fx$n); z[c(1, 3, 5, 7)] <- 1
  V_d <- randomization_cov_d(fx$X, fx$strata, fx$n1_per_stratum)
  res <- sigma_x_test(fx$X, fx$strata, z,
                      sigma_x = V_d,
                      null = "satterthwaite_asymptotic")

  ## Compute d^2 the slow way for comparison
  pi_i <- ave(z, fx$strata)
  d <- drop((z - pi_i) %*% fx$X)
  d2 <- drop(d %*% solve(V_d) %*% d)
  expect_equal(res$statistic, d2, tolerance = 1e-10)

  ## And under sigma_x = V_d, all eigenvalues are 1 so the asymptotic
  ## Satterthwaite p-value is exactly chi^2_p.
  expect_equal(res$p.value,
               pchisq(d2, df = fx$p, lower.tail = FALSE),
               tolerance = 1e-10)
})

test_that("sigma_x_test p-value is invariant under positive scaling of sigma_x", {
  ## Multiplying sigma_x by c > 0 multiplies T by 1/c and divides every
  ## lambda_k by c, so the integrand in P(sum lambda chi^2_1 > T) is unchanged.
  ## The p-value must be identical under any positive rescaling.
  fx <- make_small_fixture()
  z  <- numeric(fx$n); z[c(1, 3, 5, 7)] <- 1
  Sx <- default_sigma_x(fx$X, fx$strata)
  for (cval in c(0.1, 1, 7.3)) {
    p1 <- sigma_x_test(fx$X, fx$strata, z, sigma_x = Sx,
                       null = "satterthwaite_asymptotic")$p.value
    p2 <- sigma_x_test(fx$X, fx$strata, z, sigma_x = cval * Sx,
                       null = "satterthwaite_asymptotic")$p.value
    expect_equal(p1, p2, tolerance = 1e-10,
                 info = paste("scale =", cval))
  }
})

test_that("sigma_x_test handles a singular sigma_x via pseudoinverse", {
  ## A perfectly collinear pair of covariates produces a singular Sigma_x.
  ## The test stat should still be finite (using the Moore-Penrose pseudo-
  ## inverse) and the effective degrees of freedom should be the rank, not
  ## the column count.
  set.seed(11); n <- 20
  x1 <- rnorm(n); x3 <- rnorm(n)
  X  <- cbind(x1, x2 = 2 * x1, x3)             # rank 2, not 3
  strata <- factor(rep("--", n))
  z <- rep(c(0, 1), n / 2)
  res <- sigma_x_test(X, strata, z,
                      sigma_x = NULL,
                      null = "satterthwaite_asymptotic")
  expect_true(is.finite(res$statistic))
  ## p-value of an asymptotic 2-df chi^2 quantity
  expect_lte(length(res$lambda), 2L)
})

## --- 5.  Finite-sample moment machinery -----------------------------------
##
## The most load-bearing tests in the file: verify that the second-order moment
## machinery (ported from i113-highermoments) reproduces the EXACT mean and
## variance of T = d' Sigma_x^{-1} d under within-stratum SRSWOR, on small
## examples where we can enumerate the full randomization distribution.

test_that("sigma_x_T_moments E[T] matches exact enumeration (equal strata)", {
  fx <- make_small_fixture()
  Sx <- default_sigma_x(fx$X, fx$strata)
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  ed <- exact_randomization_dist(zs, fx$X, fx$strata, sigma_x = Sx)
  mom <- sigma_x_T_moments(fx$X, fx$strata, fx$n1_per_stratum, sigma_x = Sx)
  expect_equal(mom$EM, mean(ed$T), tolerance = 1e-12)
})

test_that("sigma_x_T_moments Var[T] matches exact enumeration (equal strata)", {
  ## This is the test that protects backend (3).  If the second-order moment
  ## formulas (Finucan / strata_t2_covariance_matrices) are off by a coefficient,
  ## this test fails.
  fx <- make_small_fixture()
  Sx <- default_sigma_x(fx$X, fx$strata)
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  ed <- exact_randomization_dist(zs, fx$X, fx$strata, sigma_x = Sx)
  exact_VarT <- mean((ed$T - mean(ed$T))^2)  # finite-population variance
  mom <- sigma_x_T_moments(fx$X, fx$strata, fx$n1_per_stratum, sigma_x = Sx)
  expect_equal(mom$VarM, exact_VarT, tolerance = 1e-10)
})

test_that("sigma_x_T_moments matches exact enumeration on uneven strata", {
  ## Hits the small-N coefficient branch (n_s = 3) and unequal n1.
  fx <- make_uneven_fixture()
  Sx <- default_sigma_x(fx$X, fx$strata)
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  ed <- exact_randomization_dist(zs, fx$X, fx$strata, sigma_x = Sx)
  mom <- sigma_x_T_moments(fx$X, fx$strata, fx$n1_per_stratum, sigma_x = Sx)
  expect_equal(mom$EM,   mean(ed$T),                       tolerance = 1e-10)
  expect_equal(mom$VarM, mean((ed$T - mean(ed$T))^2),      tolerance = 1e-10)
})

test_that("sigma_x_T_moments matches enumeration with a non-default sigma_x", {
  ## Independent of which Sigma_x the user supplies, the moment formulas must
  ## give the right E[T] and Var[T] for that specific Sigma_x.
  fx <- make_small_fixture()
  custom_Sx <- diag(c(0.5, 1, 2)) +
                 0.25 * (matrix(1, fx$p, fx$p) - diag(fx$p))
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  ed <- exact_randomization_dist(zs, fx$X, fx$strata, sigma_x = custom_Sx)
  mom <- sigma_x_T_moments(fx$X, fx$strata, fx$n1_per_stratum,
                           sigma_x = custom_Sx)
  expect_equal(mom$EM,   mean(ed$T),                  tolerance = 1e-10)
  expect_equal(mom$VarM, mean((ed$T - mean(ed$T))^2), tolerance = 1e-10)
})

## --- 6.  Backends compared: agreement and disagreement --------------------

test_that("at small n, asymptotic Var[T] OVERSTATES exact randomization Var[T]", {
  ## This pins down the substantive finding that motivated backend (3): the
  ## Gaussian identity Var[T] = 2 sum lambda^2 systematically overstates the
  ## finite-sample variance.  A future maintainer who "fixes" this collapse
  ## (e.g. by deleting backend (3) and using only the asymptotic backend) will
  ## fail this test.
  fx <- make_small_fixture()
  Sx <- default_sigma_x(fx$X, fx$strata)
  V_d <- randomization_cov_d(fx$X, fx$strata, fx$n1_per_stratum)

  ## Asymptotic Var[T] from eigenvalues of Sigma_x^{-1} V_d
  M <- solve(Sx) %*% V_d
  lambda <- Re(eigen(M, only.values = TRUE)$values)
  Var_asymp <- 2 * sum(lambda^2)

  ## Exact Var[T] from enumeration
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  ed <- exact_randomization_dist(zs, fx$X, fx$strata, sigma_x = Sx)
  Var_exact <- mean((ed$T - mean(ed$T))^2)

  expect_gt(Var_asymp, Var_exact)              # asymptotic is larger
  expect_gt(Var_asymp / Var_exact, 1.5)        # by a meaningful factor
})

test_that("draw_within_stratum_z preserves stratum counts (simple before/after)", {
  ## Simple direct check, in plain English: take an observed z, count treated
  ## and controls in each stratum, apply the permutation function, count
  ## again, and verify (a) the counts are the same and (b) the labels
  ## actually shifted (i.e. we got a real permutation, not the identity).
  set.seed(20260417)
  strata <- factor(c(rep("a", 4), rep("b", 6)))
  z_orig <- c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0)   # 2 treated in a, 3 in b

  n1_before <- tapply(z_orig,     strata, sum)
  n0_before <- tapply(1 - z_orig, strata, sum)

  by_s <- split(seq_along(z_orig), strata)
  z_new <- draw_within_stratum_z(by_s,
                                 n1_per_stratum = as.integer(n1_before),
                                 N = length(z_orig))

  n1_after <- tapply(z_new,     strata, sum)
  n0_after <- tapply(1 - z_new, strata, sum)

  expect_equal(n1_after, n1_before)
  expect_equal(n0_after, n0_before)
  expect_true(any(z_new != z_orig))           # the permutation moved units
})

test_that("draw_within_stratum_z preserves per-stratum treated counts on every draw", {
  ## Direct structural test of the simulator.  For every random draw, the
  ## number of treated units in each stratum must equal n1_per_stratum --
  ## always, no exceptions, no Monte Carlo tolerance.  This protects against
  ## a future refactor breaking the within-stratum SRSWOR property of the
  ## sampler in a way that the indirect distribution-comparison tests would
  ## not catch quickly.
  set.seed(20260414)
  ## Three strata of unequal sizes and unequal treated counts
  strata <- factor(c(rep("a", 5), rep("b", 7), rep("c", 4)))
  n1_per_stratum <- c(a = 2L, b = 3L, c = 1L)
  N <- length(strata)
  by_s <- split(seq_len(N), strata)

  ## Run many draws and check the structural invariant on every one
  B <- 500L
  for (b in seq_len(B)) {
    z <- draw_within_stratum_z(by_s, n1_per_stratum, N)
    counts <- as.integer(tapply(z, strata, sum))
    expect_equal(counts, as.integer(n1_per_stratum),
                 info = paste("draw", b))
    ## Also check z is binary 0/1
    expect_true(all(z %in% c(0, 1)))
    ## And that the total treated count matches the sum
    expect_equal(sum(z), sum(n1_per_stratum))
  }
})

test_that("draw_within_stratum_z is uniform over the SRSWOR sample space", {
  ## A 2-stratum sanity check on uniformity: every (stratum_a_assignment,
  ## stratum_b_assignment) pair should appear with frequency 1/(C(4,2)*C(4,2))
  ## = 1/36 in the limit.  We use a chi-square goodness-of-fit test against
  ## the uniform on the 36-cell joint table.
  set.seed(20260415)
  strata <- factor(rep(c("a", "b"), each = 4))
  n1_per_stratum <- c(a = 2L, b = 2L)
  N <- length(strata)
  by_s <- split(seq_len(N), strata)

  B <- 36000L
  signatures <- replicate(B, {
    z <- draw_within_stratum_z(by_s, n1_per_stratum, N)
    paste(z, collapse = "")
  })
  tab <- table(signatures)
  ## There are 6 * 6 = 36 distinct signatures; we should see all of them
  expect_equal(length(tab), 36L)
  ## Goodness of fit against uniform expectation B/36 per cell
  expected <- rep(B / 36, 36)
  observed <- as.numeric(tab)
  chisq <- sum((observed - expected)^2 / expected)
  ## Under uniform sampling, chisq ~ chi^2_{35}.  Reject if p < 0.001.
  pval <- pchisq(chisq, df = 35, lower.tail = FALSE)
  expect_gt(pval, 0.001)
})

test_that("simulate backend agrees with exact enumeration on a small example", {
  ## With B = number of enumerated assignments, the simulation backend should
  ## hit the exact distribution (when the simulator draws uniformly from the
  ## same combinatorial space).  We use n_simulate = 10000 and a moderate
  ## tolerance for the empirical CDF agreement.
  skip_on_cran()
  set.seed(20260411)
  fx <- make_small_fixture()
  Sx <- default_sigma_x(fx$X, fx$strata)
  z  <- numeric(fx$n); z[c(1, 3, 5, 7)] <- 1
  res <- sigma_x_test(fx$X, fx$strata, z,
                      sigma_x = Sx,
                      null = "simulate", n_simulate = 10000)
  ## The simulation backend should expose its draws or summaries
  expect_true(is.numeric(res$null_draws))
  expect_length(res$null_draws, 10000)

  ## Empirical mean of simulated T should match exact E[T] within MC error
  zs <- enumerate_strata_assignments(fx$strata, fx$n1_per_stratum)
  ed <- exact_randomization_dist(zs, fx$X, fx$strata, sigma_x = Sx)
  expect_equal(mean(res$null_draws), mean(ed$T),
               tolerance = 0.05 * mean(ed$T))
})

test_that("at large n, all four backends agree to within Monte Carlo error", {
  ## The asymptotic regime: when n is large, asymptotic Satterthwaite, Imhof,
  ## finite-sample Satterthwaite, and simulate must all give very close
  ## p-values.  This test confirms the asymptotic backends aren't *wrong*;
  ## they're just inappropriate at small n.
  skip_on_cran()
  skip_if_not_installed("CompQuadForm")
  set.seed(20260412)
  n <- 200; p <- 2; K <- 10
  strata <- factor(rep(seq_len(K), each = n / K))
  X <- matrix(rnorm(n * p), n, p)
  z <- unlist(lapply(seq_len(K), function(s) {
    out <- rep(0, n / K); out[sample(n / K, n / K / 2)] <- 1; out
  }))

  Sx <- default_sigma_x(X, strata)
  p_imhof <- sigma_x_test(X, strata, z, sigma_x = Sx,
                          null = "imhof")$p.value
  p_satS  <- sigma_x_test(X, strata, z, sigma_x = Sx,
                          null = "satterthwaite_asymptotic")$p.value
  p_satF  <- sigma_x_test(X, strata, z, sigma_x = Sx,
                          null = "satterthwaite_finite")$p.value
  p_sim   <- sigma_x_test(X, strata, z, sigma_x = Sx,
                          null = "simulate", n_simulate = 5000)$p.value
  expect_equal(p_imhof, p_satS, tolerance = 0.01)
  expect_equal(p_imhof, p_satF, tolerance = 0.05)
  expect_equal(p_imhof, p_sim,  tolerance = 0.05)
})

## --- 7.  End-to-end through balanceTest() ---------------------------------
##
## The new stat lives alongside the existing d^2 in the same `overall` table.
## When sigma_x_test = TRUE, three new columns ("sigma_x", "sigma_x.df",
## "sigma_x.p.value") are appended to bt$overall while the original three
## columns ("chisquare", "df", "p.value") are preserved.  This lets us read
## both stats off the same row (one row per stratification) during development.
SIGMA_X_COLS <- c("sigma_x", "sigma_x.df", "sigma_x.p.value")

test_that("balanceTest(sigma_x_test = TRUE) appends new columns to overall", {
  data(nuclearplants, package = "RItools")
  bt <- balanceTest(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
                    data = nuclearplants,
                    sigma_x_test = TRUE)
  expect_s3_class(bt, "balancetest")
  ## Original columns still present and first
  expect_true(all(c("chisquare", "df", "p.value") %in% colnames(bt$overall)))
  ## New columns appended
  expect_true(all(SIGMA_X_COLS %in% colnames(bt$overall)))
  ## Values are finite, df positive, p in [0,1]
  expect_true(all(is.finite(bt$overall[, "sigma_x"])))
  expect_true(all(bt$overall[, "sigma_x.df"] > 0))
  expect_true(all(bt$overall[, "sigma_x.p.value"] >= 0 &
                  bt$overall[, "sigma_x.p.value"] <= 1))
})

test_that("balanceTest(sigma_x_test = FALSE) is byte-identical to current behavior", {
  ## Backwards compatibility: opt-in must mean opt-in.  Without the new flag
  ## the result has the historical columns and only the historical columns.
  data(nuclearplants, package = "RItools")
  bt_off <- balanceTest(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
                        data = nuclearplants)
  expect_named(bt_off$overall, c("chisquare", "df", "p.value"),
               ignore.order = FALSE)
  expect_false(any(SIGMA_X_COLS %in% colnames(bt_off$overall)))
})

test_that("balanceTest existing chisquare column is unchanged when sigma_x_test = TRUE", {
  ## Adding the new columns must NOT alter the existing d^2 numbers.  This
  ## catches accidental refactors that recompute HB08 in a different way.
  data(nuclearplants, package = "RItools")
  fmla <- pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n
  bt_off <- balanceTest(fmla, data = nuclearplants)
  bt_on  <- balanceTest(fmla, data = nuclearplants, sigma_x_test = TRUE)
  expect_equal(bt_on$overall[, c("chisquare", "df", "p.value")],
               bt_off$overall[, c("chisquare", "df", "p.value")],
               tolerance = 1e-12)
})

test_that("balanceTest sigma_x stat agrees with manual sigma_x_test() call", {
  ## End-to-end correctness check: bt$overall["pt", "sigma_x"] equals the
  ## number you would get by extracting d, strata, and Sigma_x by hand and
  ## calling the mid-level sigma_x_test() function directly.
  data(nuclearplants, package = "RItools")
  bt <- balanceTest(pr ~ date + t1 + t2 + strata(pt),
                    data = nuclearplants,
                    sigma_x_test = TRUE,
                    null = "satterthwaite_finite")

  X <- as.matrix(nuclearplants[, c("date", "t1", "t2")])
  strata <- factor(nuclearplants$pt)
  z <- nuclearplants$pr
  Sx <- default_sigma_x(X, strata)
  manual <- sigma_x_test(X, strata, z, sigma_x = Sx,
                         null = "satterthwaite_finite")
  expect_equal(bt$overall["pt", "sigma_x"],
               manual$statistic, tolerance = 1e-8)
  expect_equal(bt$overall["pt", "sigma_x.p.value"],
               manual$p.value,   tolerance = 1e-8)
})

test_that("balanceTest 'null' argument is honored", {
  ## The user must be able to pick the null backend at the top level.  Two
  ## different choices ('satterthwaite_finite' and 'satterthwaite_asymptotic')
  ## should give DIFFERENT p-values on a small enough design (where the
  ## asymptotic-vs-finite gap is real); they would coincide as n grows.
  data(nuclearplants, package = "RItools")
  bt_finite <- balanceTest(pr ~ date + t1 + t2 + strata(pt),
                           data = nuclearplants, sigma_x_test = TRUE,
                           null = "satterthwaite_finite")
  bt_asymp  <- balanceTest(pr ~ date + t1 + t2 + strata(pt),
                           data = nuclearplants, sigma_x_test = TRUE,
                           null = "satterthwaite_asymptotic")
  ## Same statistic, different p-values
  expect_equal(bt_finite$overall[, "sigma_x"],
               bt_asymp$overall[, "sigma_x"], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(bt_finite$overall[, "sigma_x.p.value"],
                                bt_asymp$overall[, "sigma_x.p.value"],
                                tolerance = 1e-4)))
})

test_that("subset.xbal preserves sigma_x columns and attribute", {
  ## Subsetting must keep the new columns and the sigma_x_info attribute, so
  ## downstream tooling that operates on a subset still has access to lambda
  ## and the sigma_x matrix per stratification.
  data(nuclearplants, package = "RItools")
  bt <- balanceTest(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
                    data = nuclearplants, sigma_x_test = TRUE)
  bt_pt <- subset(bt, strata = "pt")
  expect_true(all(SIGMA_X_COLS %in% colnames(bt_pt$overall)))
  expect_true(!is.null(attr(bt_pt$overall, "sigma_x_info")))
  expect_equal(names(attr(bt_pt$overall, "sigma_x_info")), "pt")
})

test_that("glance.xbal includes sigma_x columns when present", {
  ## Confirm the broom-style glance method picks up the new columns with no
  ## additional code.  This is a smoke test for the column-pass-through claim
  ## in the glance.xbal docstring.
  data(nuclearplants, package = "RItools")
  bt <- balanceTest(pr ~ date + t1 + t2 + strata(pt),
                    data = nuclearplants, sigma_x_test = TRUE)
  g <- glance.xbal(bt, strata = "pt")
  expect_true(all(SIGMA_X_COLS %in% colnames(g)))
  expect_equal(nrow(g), 1L)
})

test_that("balanceTest sigma_x stat flows through clusters (smoke)", {
  ## With a clustered design, the new test stat should at least run end-to-
  ## end without error and return finite values.  Exact numerical agreement
  ## with cluster-aggregated manual computation is covered by the moment-
  ## machinery tests on raw inputs above.
  data(ym_long, package = "RItools")
  bt <- balanceTest(trt ~ assessed + hypo + lipid +
                          strata(assess_strata) + cluster(practice),
                    data = ym_long,
                    sigma_x_test = TRUE)
  expect_true(all(is.finite(bt$overall[, "sigma_x"])))
  expect_true(all(bt$overall[, "sigma_x.p.value"] >= 0 &
                  bt$overall[, "sigma_x.p.value"] <= 1))
})
