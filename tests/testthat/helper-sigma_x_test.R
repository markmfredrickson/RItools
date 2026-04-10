## Helpers for tests of the sigma_x balance test (d' Sigma_x^{-1} d).
##
## These give exact-randomization ground truth on tiny examples by enumerating
## every assignment compatible with the within-stratum treated counts.  The
## moment- and CDF-based backends in the package are then verified against
## these exact answers.  Adapted from misc.R on the i113-highermoments branch
## (Mark Fredrickson).
##
## Use only for *small* designs.  The number of enumerated assignments is
## prod_s choose(n_s, n1_s); blows up fast.

## Enumerate every treatment assignment compatible with within-stratum SRSWOR.
##
## strata          factor of length n
## n1_per_stratum  integer vector of length nlevels(strata)
##
## Returns an n x B matrix; column b is one z vector.
enumerate_strata_assignments <- function(strata, n1_per_stratum) {
  stopifnot(is.factor(strata),
            length(n1_per_stratum) == nlevels(strata))
  n <- length(strata)
  K <- nlevels(strata)
  per_s <- vector("list", K)
  for (s in seq_len(K)) {
    idx_in_s <- which(as.integer(strata) == s)
    n_s <- length(idx_in_s)
    if (n1_per_stratum[s] < 1L || n1_per_stratum[s] >= n_s) {
      stop("each stratum needs at least 1 treated and 1 control")
    }
    per_s[[s]] <- list(idx = idx_in_s,
                       combos = combn(n_s, n1_per_stratum[s]))
  }
  ncol_per_s <- vapply(per_s, function(p) ncol(p$combos), integer(1))
  combo_idx  <- as.matrix(do.call(expand.grid, lapply(ncol_per_s, seq_len)))
  B <- nrow(combo_idx)
  out <- matrix(0, nrow = n, ncol = B)
  for (b in seq_len(B)) {
    for (s in seq_len(K)) {
      ps <- per_s[[s]]
      treated_locs <- ps$combos[, combo_idx[b, s]]
      out[ps$idx[treated_locs], b] <- 1
    }
  }
  out
}

## Exact randomization distribution of d_raw and T = d' Sigma_x^{-1} d under
## within-stratum SRSWOR.  Returns the per-assignment d vectors, the per-
## assignment T values, and the empirical Cov(d) (finite-population variance).
##
## d_raw is computed as (z - pi)' X with pi the within-stratum treated
## proportion.  Because (z - pi) is mean zero within each stratum,
## (z - pi)' (X - S Xbar) = (z - pi)' X, so the d we report is independent of
## any stratum-mean centering.
exact_randomization_dist <- function(zs, X, strata, sigma_x) {
  stopifnot(is.matrix(zs),
            nrow(zs) == nrow(X),
            nrow(zs) == length(strata))
  n_per_s  <- tabulate(as.integer(strata))
  n1_per_s <- tapply(zs[, 1], strata, sum)
  pi_i     <- (n1_per_s / n_per_s)[as.integer(strata)]

  ## d = (z - pi)' X for every column of zs at once: p x B
  zc    <- sweep(zs, 1, pi_i, "-")
  d_mat <- t(crossprod(zc, X))

  Si    <- solve(sigma_x)
  T_vec <- vapply(seq_len(ncol(d_mat)), function(b) {
    drop(d_mat[, b] %*% Si %*% d_mat[, b])
  }, numeric(1))

  ## Use 1/B (finite-population) rather than 1/(B-1) so that the empirical
  ## moments equal the exact randomization moments to machine precision.
  B <- ncol(d_mat)
  V_d_emp <- crossprod(t(d_mat)) / B - tcrossprod(rowMeans(d_mat))

  list(d = d_mat, T = T_vec, V_d_emp = V_d_emp)
}

## Closed form for Cov(d_raw) under within-stratum SRSWOR:
##   Cov(d_raw) = sum_s [ n1_s n0_s / (n_s (n_s - 1)) ] * S_xs
## where S_xs = sum_{i in s} (x_i - xbar_s)(x_i - xbar_s)' is the unscaled
## within-stratum sum of squared deviations.  Verified empirically in
## test_that("randomization_cov_d closed form ...", ...).
closed_form_V_d_helper <- function(X, strata, n1_per_stratum) {
  K <- nlevels(strata)
  stopifnot(length(n1_per_stratum) == K)
  n_per_s <- tabulate(as.integer(strata))
  by_s    <- split(seq_len(nrow(X)), strata)
  V <- matrix(0, ncol(X), ncol(X))
  for (s in seq_len(K)) {
    idx <- by_s[[s]]
    n_s <- n_per_s[s]
    n1_s <- n1_per_stratum[s]
    n0_s <- n_s - n1_s
    if (n_s < 2) next
    Xc <- scale(X[idx, , drop = FALSE], scale = FALSE)
    S_xs <- crossprod(Xc)
    V <- V + (n1_s * n0_s / (n_s * (n_s - 1))) * S_xs
  }
  V
}

## Default within-stratum-pooled sample covariance, the proposed default
## Sigma_x:
##   Sigma_x_default = (1/(N - K)) * sum_s sum_{i in s} (x_i - xbar_s)(x_i - xbar_s)'
## This is the unbiased estimator from a stratum-fixed-effects ANOVA residual,
## and it is the multivariate analogue of the existing pooled-SD denominator
## that balanceTest uses for the univariate std.diff column.
within_stratum_pooled_cov_helper <- function(X, strata) {
  K <- nlevels(strata)
  N <- nrow(X)
  by_s <- split(seq_len(N), strata)
  S <- matrix(0, ncol(X), ncol(X))
  for (idx in by_s) {
    if (length(idx) < 2) next
    Xc <- scale(X[idx, , drop = FALSE], scale = FALSE)
    S <- S + crossprod(Xc)
  }
  S / (N - K)
}

## A small canonical test fixture: 2 strata of size 4, 2 treated per stratum,
## 3 covariates, set seed for reproducibility.  6 * 6 = 36 enumerated z's.
make_small_fixture <- function(seed = 20260409) {
  set.seed(seed)
  n <- 8L; K <- 2L; p <- 3L
  list(n = n, K = K, p = p,
       strata = factor(rep(seq_len(K), each = n / K)),
       n1_per_stratum = c(2L, 2L),
       X = matrix(rnorm(n * p), n, p,
                  dimnames = list(NULL, paste0("x", seq_len(p)))))
}

## A second fixture with unequal stratum sizes (3 + 5 = 8) and unequal treated
## counts (1, 2).  10 * 10 = 100 enumerated z's.  Hits the small-N coefficient
## branches in the moment machinery.
make_uneven_fixture <- function(seed = 20260410) {
  set.seed(seed)
  n <- 8L; p <- 2L
  list(n = n, K = 2L, p = p,
       strata = factor(c(rep(1L, 3L), rep(2L, 5L))),
       n1_per_stratum = c(1L, 2L),
       X = matrix(rnorm(n * p), n, p,
                  dimnames = list(NULL, paste0("x", seq_len(p)))))
}
