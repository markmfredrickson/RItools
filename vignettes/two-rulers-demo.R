## ---------------------------------------------------------------------------
## two-rulers-demo.R
##
## Jake's exact proposal, and exactly where it breaks.
##
## Proposal: don't pick 0.25 out of the air. Set the max|SMD| yardstick from the
## max|SMD| you would see across 1000 re-randomizations of the ANALOGOUS BLOCK
## RANDOMIZED EXPERIMENT (re-randomize treatment within the matched sets).
##
## The problem: there are TWO rulers for "how big is the imbalance," and they
## disagree by construction as the match tightens.
##   Ruler A (FIXED): the covariate's own units, or its fixed pooled SD. A 0.2-unit
##     gap is 0.2 units whether the strata are loose or tight. max|SMD| with a
##     fixed pooled-SD denominator is Ruler A. It does NOT move as you tighten.
##   Ruler B (SHRINKS): the spread of max|SMD| under within-set re-randomization,
##     which is set by the within-set randomization variance V_d. Tight sets leave
##     z almost no room to move imbalance, so this spread -> 0 as the match tightens.
## Calibrating max|SMD| against within-set re-randomization MEASURES THE FIXED GAP
## WITH RULER B. The reference shrinks under the observed value, the percentile
## climbs to 1, and the omnibus p falls -- even though the gap in covariate units
## never changed. The escape (Ruler A): a FIXED reference -- a substantively chosen
## tolerance, or re-randomization within a FIXED, COARSER design that does not
## shrink with the match.
##
## We hold a persistent per-covariate gap g FIXED (raw, in covariate units) and
## TIGHTEN the match (shrink within-set noise). Fixed units would do the same; the
## shrink knob just makes the collapse clean and isolates V_d.
## ---------------------------------------------------------------------------

set.seed(20260617)

## ---- machinery (matches randomization_cov_d / balanceTestEngine) ----------
adj_diff_sum <- function(X, strata, z) {     # d = (z - pi)' X
  pi <- ave(z, strata, FUN = mean)
  drop(crossprod(z - pi, as.matrix(X)))
}
Vperm <- function(X, strata, z) {            # within-block randomization cov of d
  X <- as.matrix(X); V <- matrix(0, ncol(X), ncol(X))
  for (idx in split(seq_len(nrow(X)), strata)) {
    n <- length(idx); n1 <- sum(z[idx]); n0 <- n - n1
    if (n < 2 || n1 == 0 || n0 == 0) next
    Xc <- scale(X[idx, , drop = FALSE], scale = FALSE)
    V  <- V + (n1 * n0 / (n * (n - 1))) * crossprod(Xc)
  }
  V
}
d2_stat <- function(d, V, tol = sqrt(.Machine$double.eps)) {
  V <- (V + t(V)) / 2; e <- eigen(V, symmetric = TRUE)
  pos <- e$values > tol * max(e$values, 1)
  if (!any(pos)) return(0)
  proj <- crossprod(e$vectors[, pos, drop = FALSE], d)
  sum(proj^2 / e$values[pos])
}
perm_within <- function(strata, z) {
  out <- z; for (idx in split(seq_along(z), strata)) out[idx] <- sample(z[idx]); out
}

## weighted within-stratum standardized mean difference vector (Ruler A: fixed
## pooled-SD denominator, computed once on the observed data). Harmonic weights.
smd_vec <- function(X, strata, z, sd_pool) {
  X <- as.matrix(X); K <- ncol(X)
  num <- numeric(K); wsum <- 0
  for (idx in split(seq_len(nrow(X)), strata)) {
    zz <- z[idx]; n1 <- sum(zz); n0 <- length(zz) - n1
    if (n1 == 0 || n0 == 0) next
    w  <- n1 * n0 / (n1 + n0)
    mt <- colMeans(X[idx, , drop = FALSE][zz == 1, , drop = FALSE])
    mc <- colMeans(X[idx, , drop = FALSE][zz == 0, , drop = FALSE])
    num <- num + w * (mt - mc); wsum <- wsum + w
  }
  (num / wsum) / sd_pool
}
pooled_sd <- function(X, z) {                # fixed pooled within-group SD (Ruler A scale)
  X <- as.matrix(X)
  apply(X, 2, function(x) {
    v1 <- var(x[z == 1]); v0 <- var(x[z == 0])
    n1 <- sum(z == 1); n0 <- sum(z == 0)
    sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
  })
}

## hierarchical data: 4 coarse groups x 4 matched sets x (2T+2C). Persistent raw
## per-covariate gap g (FIXED). Within-set noise scaled by `shrink` (the match).
make_data <- function(shrink, g = 0.2, K = 2,
                      group_centers = c(0, 10, 20, 30), sets_per_group = 4) {
  Xlist <- list(); ms <- integer(0); grp <- integer(0); z <- integer(0); setid <- 0L
  for (gci in seq_along(group_centers)) for (s in seq_len(sets_per_group)) {
    setid <- setid + 1L
    ctr <- group_centers[gci] + rnorm(K, sd = 1.5)
    Xlist[[setid]] <- rbind(ctr + g/2 + shrink*rnorm(K), ctr + g/2 + shrink*rnorm(K),
                            ctr - g/2 + shrink*rnorm(K), ctr - g/2 + shrink*rnorm(K))
    ms <- c(ms, rep(setid, 4)); grp <- c(grp, rep(gci, 4)); z <- c(z, c(1,1,0,0))
  }
  X <- do.call(rbind, Xlist); colnames(X) <- paste0("x", seq_len(K))
  list(X = X, set = factor(ms), group = factor(grp), z = z)
}

## re-randomization distribution of max|SMD| under a given reference blocking
rerand_maxsmd <- function(X, ref_strata, z, sd_pool, B = 2000) {
  obs <- max(abs(smd_vec(X, ref_strata, z, sd_pool)))
  null <- numeric(B)
  for (b in seq_len(B)) {
    zp <- perm_within(ref_strata, z)
    null[b] <- max(abs(smd_vec(X, ref_strata, zp, sd_pool)))
  }
  list(obs = obs, q95 = quantile(null, 0.95),
       pctile = mean(null <= obs - 1e-12))         # fraction of refs at or below obs
}

## d^2 omnibus p at a blocking (analytic chisq, HB08)
d2_p <- function(X, strata, z) {
  d <- adj_diff_sum(X, strata, z); V <- Vperm(X, strata, z)
  stat <- d2_stat(d, V)
  rank <- sum(eigen((V + t(V))/2, symmetric = TRUE)$values >
                sqrt(.Machine$double.eps) * max(eigen((V+t(V))/2, symmetric=TRUE)$values, 1))
  pchisq(stat, df = rank, lower.tail = FALSE)
}

## ---- run across tightening matches ----------------------------------------
shrinks <- c(0.60, 0.30, 0.10, 0.03, 0.01)
tab <- data.frame()
for (i in seq_along(shrinks)) {
  set.seed(400 + i)
  d <- make_data(shrinks[i])
  sdp <- pooled_sd(d$X, d$z)
  ## raw gap in covariate units (Ruler A, no SD at all): the weighted within-set
  ## mean difference, max over covariates, in the covariate's own units
  raw <- max(abs(smd_vec(d$X, d$set, d$z, sd_pool = rep(1, ncol(d$X)))))
  ## Ruler B: re-randomize within the MATCHED SETS (the "analogous block experiment")
  setref   <- rerand_maxsmd(d$X, d$set,   d$z, sdp)
  ## A FIXED, coarser reference: re-randomize within the 4 COARSE GROUPS
  groupref <- rerand_maxsmd(d$X, d$group, d$z, sdp)
  tab <- rbind(tab, data.frame(
    shrink          = shrinks[i],
    raw_max_gap     = raw,                       # Ruler A, covariate units: ~ fixed
    obs_max_smd     = setref$obs,                # Ruler A, pooled-SD: ~ fixed
    set_rerand_q95  = as.numeric(setref$q95),    # Ruler B threshold: SHRINKS
    set_pctile_obs  = setref$pctile,             # -> 1 (collapse)
    d2_p_set        = d2_p(d$X, d$set, d$z),      # -> 0 (collapse)
    group_rerand_q95 = as.numeric(groupref$q95), # FIXED reference: stable
    group_pctile_obs = groupref$pctile           # stable verdict
  ))
}

cat("==== Two rulers for the SAME fixed imbalance, as the match tightens ====\n\n")
cat("Ruler A (fixed): raw_max_gap (covariate units) and obs_max_smd (pooled SD).\n")
cat("Ruler B (shrinks): set_rerand_q95 = 95th pct of max|SMD| under within-SET\n")
cat("   re-randomization -- the yardstick Jake proposed. set_pctile_obs and\n")
cat("   d2_p_set are the observed read against it.\n")
cat("FIXED coarse reference: group_rerand_q95 / group_pctile_obs (within-GROUP).\n\n")
op <- options(digits = 3); print(tab, row.names = FALSE); options(op)

cat("\nReading: raw_max_gap and obs_max_smd barely move (the imbalance is FIXED).\n")
cat("set_rerand_q95 collapses toward 0, so set_pctile_obs climbs toward 1 and\n")
cat("d2_p_set toward 0 -- Jake's re-randomized max|SMD| yardstick FAILS for the\n")
cat("same reason the omnibus p does: it is Ruler B. The within-GROUP reference is\n")
cat("FIXED and gives a stable verdict -- Ruler A -- but it compares to a coarser,\n")
cat("UNBUILT design (the trilemma cost).\n")

saveRDS(tab, "vignettes/two-rulers-demo.rds")
cat("\ndata written to vignettes/two-rulers-demo.rds\n")
