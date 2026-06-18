## Shared helpers for the balance-collapse / screening / ACAT simulation studies.
## Sourced by sim-collapse-pvalue-ess.R and sim-screening-acat.R.
##
## Design point (see fixed-metric-collapse-memo.md and
## screening-and-pvalue-ess-memo.md): the within-set d^2 omnibus is homogeneous
## of degree zero in the within-pair differences, so it reads the DIRECTION/sign
## of imbalance relative to within-pair reshuffling and is blind to MAGNITUDE.
## These helpers generate matched-pair data in which we can dial, separately,
## (i) the systematic directional bias mu (the alternative the test has power
## against), (ii) the within-pair dissimilarity sigma_c (which sets V_d, i.e. the
## "magnitude" the test is blind to), (iii) the number of pairs S (= effective
## sample size for pairs), and (iv) a global scale t applied to BOTH mu and
## sigma_c (the pure degree-zero / collapse knob).

## Load the package whether this file is sourced from the repo root (the .R
## driver scripts: DESCRIPTION is in ".") or from vignettes/ (the .qmd: the
## package root is "..").
.pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."
suppressMessages(devtools::load_all(.pkg_root, export_all = FALSE, quiet = TRUE))

## ---- data-generating process ------------------------------------------------
## For each pair s the two units share a pair location m_s (between-pair
## variation -> what stratification removes) and differ by
##     treated - control = D_s = noise_s + mu,   noise_s ~ N(0, sigma_c^2),
## with the treated unit being the (+D_s/2) side.  mu is the SYSTEMATIC part
## (constant across pairs, in expectation the only thing dbar sees); noise_s is
## symmetric, so under the within-pair permutation null its sign is exchangeable.
## Multiplying sigma_c and mu by `tscale` shrinks the whole imbalance without
## changing its direction -- the degree-zero knob.
##
## sigma_c and mu are length-k vectors.  An exactly/near-exactly matched
## covariate is just one with sigma_c[j] ~ 0 and mu[j] = 0.
make_pairs <- function(S, sigma_c, mu, tscale = 1, rho = 0.5) {
  k <- length(sigma_c)
  stopifnot(length(mu) == k)
  ## Covariates are CORRELATED across columns (not independent): moderate
  ## equicorrelation R[j,l] = rho off the diagonal, 1 on it (rho = 0 would
  ## recover independent covariates).  We apply the same correlation to the
  ## between-pair locations and to the within-pair differences; it is the
  ## correlation of the within-pair DIFFERENCES that puts off-diagonal entries
  ## into V_d and so actually reaches the tests.  Column scaling by sigma_c keeps
  ## each covariate's marginal within-pair SD equal to sigma_c[j], so the
  ## mu/sigma_c ("standardized bias") interpretation is unchanged.
  R <- matrix(rho, k, k); diag(R) <- 1
  L <- t(chol(R))                                              # L %*% z has cov R
  cor_rows <- function() t(L %*% matrix(rnorm(k * S), nrow = k))  # S x k, row cov R
  m     <- cor_rows()                                          # pair locations
  noise <- cor_rows() *
           matrix(sigma_c * tscale, nrow = S, ncol = k, byrow = TRUE)
  D     <- noise + matrix(mu * tscale, nrow = S, ncol = k, byrow = TRUE)  # tr - ctl
  treated <- m + D / 2
  control <- m - D / 2
  X <- matrix(0, nrow = 2 * S, ncol = k)
  X[seq(1, 2 * S, by = 2), ] <- treated                        # odd rows treated
  X[seq(2, 2 * S, by = 2), ] <- control                        # even rows control
  out <- data.frame(z = rep(c(1L, 0L), times = S),
                    pair = factor(rep(seq_len(S), each = 2)))
  for (j in seq_len(k)) out[[paste0("x", j)]] <- X[, j]
  out
}

covariate_names <- function(df) grep("^x[0-9]+$", names(df), value = TRUE)

bt_formula <- function(covs) {
  reformulate(c(covs, "strata(pair)"), response = "z")
}

## ---- run the omnibus/omnibuses on one data set ------------------------------
## Returns the d^2 chisquare, its df, the d^2 p-value, and the ACAT p-value for
## the "pair" stratification.  Wrapped in tryCatch so the degenerate-covariance
## stop() (R/utils.R:381) becomes a recorded NA + error flag instead of halting.
run_omnibus <- function(df, covs = covariate_names(df)) {
  if (length(covs) == 0L) {
    return(list(chisq = NA_real_, df = 0L, p_d2 = NA_real_,
                p_acat = NA_real_, errored = FALSE, abstained = TRUE))
  }
  res <- tryCatch({
    bt <- balanceTest(bt_formula(covs), data = df, cauchy.combination = TRUE)
    ov <- bt$overall
    list(chisq = ov["pair", "chisquare"], df = ov["pair", "df"],
         p_d2 = ov["pair", "p.value"], p_acat = ov["pair", "cauchy_comb_p"],
         errored = FALSE, abstained = FALSE)
  }, error = function(e) {
    list(chisq = NA_real_, df = NA_integer_, p_d2 = NA_real_,
         p_acat = NA_real_, errored = TRUE, abstained = FALSE)
  })
  res
}

## ---- design-based screening statistic ---------------------------------------
## w_j = within-stratum variance / total variance of covariate j: the fraction
## of the covariate's variance that survived stratification.  Exact matching
## drives it to ~0.  This is a function of the covariate values and strata ONLY
## -- NOT of the treatment vector z -- so screening on it cannot be gamed by
## peeking at the observed imbalance.  A production screen should read
## diag(tcov) / complete-rand variance from the engine; w_j is the base-R proxy.
within_var_fraction <- function(df, covs = covariate_names(df)) {
  vapply(covs, function(cn) {
    x <- df[[cn]]
    within_resid <- x - ave(x, df$pair)                # deviation from pair mean
    within_var <- sum(within_resid^2) / (length(x) - nlevels(df$pair))
    total_var  <- var(x)
    if (total_var <= 0) 0 else within_var / total_var
  }, numeric(1))
}

## Survivors of the screen at tolerance tau (keep covariates that retain at least
## a fraction tau of their variance within strata).
screen_survivors <- function(df, tau, covs = covariate_names(df)) {
  w <- within_var_fraction(df, covs)
  covs[w >= tau]
}

## ---- direction (singular-value) screen --------------------------------------
## Generalize w_j from single covariates to DIRECTIONS.  The fraction of variance
## surviving stratification in a direction v is (v' W v)/(v' Tot v), with W the
## pooled within-stratum covariance of X and Tot its total covariance.  The
## generalized eigenpairs of (W, Tot) give one fraction per orthogonal direction;
## an exactly-matched direction has fraction ~0.  Under correlated/collinear
## covariates that degenerate direction can be a LINEAR COMBINATION that no single
## covariate's w_j reveals (e.g. x1 and x2 each vary within pairs but their
## difference is matched).  We keep directions with fraction >= tau and return the
## covariates PROJECTED onto them, which balanceTest can then test directly --
## projection is linear, so the engine's tcov on the projection equals the
## restriction of the original tcov to the kept subspace.
sv_screen_project <- function(df, tau, covs = covariate_names(df)) {
  X    <- as.matrix(df[, covs, drop = FALSE])
  pair <- df$pair
  resid <- X - apply(X, 2, function(col) ave(col, pair))        # deviations from pair means
  W   <- crossprod(resid) / (nrow(X) - nlevels(pair))           # pooled within-stratum cov
  Tot <- stats::cov(X)                                          # total cov
  ## whiten by Tot, eigendecompose the symmetric M = Tot^{-1/2} W Tot^{-1/2}
  e_tot <- eigen(Tot, symmetric = TRUE)
  Tot_ih <- e_tot$vectors %*% diag(1 / sqrt(pmax(e_tot$values, 1e-12)), length(e_tot$values)) %*%
            t(e_tot$vectors)
  e_m    <- eigen(Tot_ih %*% W %*% Tot_ih, symmetric = TRUE)
  lambda <- pmax(e_m$values, 0)                                 # within-stratum fraction per direction
  keep   <- which(lambda >= tau)
  V_keep <- Tot_ih %*% e_m$vectors[, keep, drop = FALSE]        # generalized eigenvectors (kept)
  list(lambda = sort(lambda, decreasing = TRUE), n_keep = length(keep),
       P = X %*% V_keep)                                        # projected covariates, one col per kept direction
}

## Run the d^2 omnibus on the direction-screen survivors.
run_sv_screen <- function(df, tau, covs = covariate_names(df)) {
  sv <- sv_screen_project(df, tau, covs)
  if (sv$n_keep == 0L)
    return(list(p_d2 = NA_real_, df = 0L, n_keep = 0L, lambda = sv$lambda))
  dproj <- data.frame(z = df$z, pair = df$pair)
  pcols <- paste0("d", seq_len(ncol(sv$P)))
  for (j in seq_along(pcols)) dproj[[pcols[j]]] <- sv$P[, j]
  r <- run_omnibus(dproj, covs = pcols)
  list(p_d2 = r$p_d2, df = r$df, n_keep = sv$n_keep, lambda = sv$lambda)
}

## ---- effective sample size (matched-set harmonic weights) -------------------
## ESS = sum_s harmonic_mean(n1_s, n0_s) = sum_s 2 n1 n0 / (n1 + n0).
## For pairs every term is 1, so ESS = number of pairs.
effective_sample_size <- function(df) {
  tab <- table(df$pair, df$z)
  n1 <- tab[, "1"]; n0 <- tab[, "0"]
  sum(2 * n1 * n0 / (n1 + n0))
}

## ---- small utilities --------------------------------------------------------
rejection_rate <- function(p, alpha = 0.05) mean(p <= alpha, na.rm = TRUE)
