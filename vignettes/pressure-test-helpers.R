## ---------------------------------------------------------------------------
## pressure-test-helpers.R   (shared harness for the impossibility pressure-test)
##
## Question: is there a DESIGN-INTERNAL, SINGLE-NUMBER balance calibration that
## does NOT collapse as the match tightens with a FIXED substantive imbalance?
##
## A "calibration" returns an ALARM score (higher = "this design looks more
## imbalanced") by comparing an observed imbalance to a REFERENCE. The reference
## is either:
##   DESIGN-INTERNAL  -- generated from the realized strata (within-set
##                       re-randomization; any function of the within-stratum V_d).
##   EXTERNAL/FIXED   -- generated from a fixed coarser design that does NOT
##                       shrink with the match (whole-pool complete randomization,
##                       fixed-Sigma Mahalanobis, a substantive tolerance).
##
## The collapse test holds the imbalance FIXED (gap g in covariate units) and
## tightens the match (shrink -> 0, so within-set V_d -> 0). A calibration
## COLLAPSES if its alarm climbs toward maximal purely from tightening -- i.e.,
## it reports a beautifully matched design as ever more imbalanced for a gap that
## never changed.
##
## Every candidate is calib(dat) -> list(alarm = <numeric, higher=worse>,
## ref = "internal"|"external"|"none"). The harness sweeps it and judges.
## ---------------------------------------------------------------------------

## ---- low-level (matches randomization_cov_d / balanceTestEngine) -----------
adj_diff_sum <- function(X, strata, z) { pi <- ave(z, strata, FUN = mean); drop(crossprod(z - pi, as.matrix(X))) }
Vperm <- function(X, strata, z) {
  X <- as.matrix(X); V <- matrix(0, ncol(X), ncol(X))
  for (idx in split(seq_len(nrow(X)), strata)) {
    n <- length(idx); n1 <- sum(z[idx]); n0 <- n - n1
    if (n < 2 || n1 == 0 || n0 == 0) next
    Xc <- scale(X[idx, , drop = FALSE], scale = FALSE)
    V  <- V + (n1 * n0 / (n * (n - 1))) * crossprod(Xc)
  }
  V
}
eig_rank <- function(V, tol = sqrt(.Machine$double.eps)) {
  V <- (V + t(V))/2; ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  sum(ev > tol * max(ev, 1))
}
d2_stat <- function(d, V, tol = sqrt(.Machine$double.eps)) {
  V <- (V + t(V))/2; e <- eigen(V, symmetric = TRUE)
  pos <- e$values > tol * max(e$values, 1); if (!any(pos)) return(0)
  proj <- crossprod(e$vectors[, pos, drop = FALSE], d); sum(proj^2 / e$values[pos])
}
perm_within <- function(strata, z) { out <- z; for (idx in split(seq_along(z), strata)) out[idx] <- sample(z[idx]); out }
pooled_sd <- function(X, z) {
  X <- as.matrix(X)
  apply(X, 2, function(x) { v1 <- var(x[z==1]); v0 <- var(x[z==0]); n1 <- sum(z==1); n0 <- sum(z==0)
                            sqrt(((n1-1)*v1 + (n0-1)*v0)/(n1+n0-2)) })
}
## weighted within-stratum standardized mean-difference vector (fixed-SD denom)
smd_vec <- function(X, strata, z, sd_pool) {
  X <- as.matrix(X); K <- ncol(X); num <- numeric(K); wsum <- 0
  for (idx in split(seq_len(nrow(X)), strata)) {
    zz <- z[idx]; n1 <- sum(zz); n0 <- length(zz) - n1; if (n1==0 || n0==0) next
    w <- n1*n0/(n1+n0)
    mt <- colMeans(X[idx,,drop=FALSE][zz==1,,drop=FALSE]); mc <- colMeans(X[idx,,drop=FALSE][zz==0,,drop=FALSE])
    num <- num + w*(mt - mc); wsum <- wsum + w
  }
  (num/wsum)/sd_pool
}

## ---- DGP: hierarchical matched data, FIXED gap, match tightness = shrink ----
## groups (coarse, fixed structure) x sets (matched) x (2T+2C). The per-covariate
## treated-minus-control gap g is FIXED in covariate units; `shrink` scales the
## within-set noise -> as shrink->0 the match becomes exact and within-set V_d->0.
make_matched <- function(shrink, g = 0.2, K = 3, rho = 0.5,
                         group_centers = c(0, 8, 16, 24), sets_per_group = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  R <- matrix(rho, K, K); diag(R) <- 1; L <- t(chol(R))
  cn <- function(n) t(L %*% matrix(rnorm(K*n), nrow = K))
  u <- rep(1, K)/sqrt(K)
  Xl <- list(); ms <- integer(0); grp <- integer(0); z <- integer(0); sid <- 0L
  for (gc in seq_along(group_centers)) for (s in seq_len(sets_per_group)) {
    sid <- sid + 1L
    ctr <- group_centers[gc] + as.numeric(cn(1))*1.5
    blk <- rbind(ctr + (g/2)*u + shrink*as.numeric(cn(1)),
                 ctr + (g/2)*u + shrink*as.numeric(cn(1)),
                 ctr - (g/2)*u + shrink*as.numeric(cn(1)),
                 ctr - (g/2)*u + shrink*as.numeric(cn(1)))
    Xl[[sid]] <- blk; ms <- c(ms, rep(sid,4)); grp <- c(grp, rep(gc,4)); z <- c(z, c(1,1,0,0))
  }
  X <- do.call(rbind, Xl); colnames(X) <- paste0("x", seq_len(K))
  list(X = X, z = z,
       set = factor(ms), group = factor(grp), pool = factor(rep(1L, length(z))),
       sd_pool = pooled_sd(X, z))
}

## ---- reference statistics & mechanisms ------------------------------------
## complete-randomization assignment on the whole pool, preserving n1/n0
perm_pool <- function(z) sample(z)
## value of a statistic T(z) under a reference mechanism -> upper-tail percentile
ref_percentile <- function(Tobs, Tfun, gen, B = 1500) {
  null <- replicate(B, Tfun(gen())); mean(null <= Tobs - 1e-12)
}

## ============================================================================
## NAMED CONTENDERS  (alarm in [0,1] unless noted; higher = looks more imbalanced)
## ============================================================================

## 1. d^2 omnibus at the matched sets (analytic chisq). DESIGN-INTERNAL. baseline.
cal_d2_internal <- function(dat) {
  d <- adj_diff_sum(dat$X, dat$set, dat$z); V <- Vperm(dat$X, dat$set, dat$z)
  p <- pchisq(d2_stat(d, V), df = max(eig_rank(V),1), lower.tail = FALSE)
  list(alarm = 1 - p, ref = "internal")
}

## 2. max|SMD| percentile vs within-SET re-randomization. DESIGN-INTERNAL. (Jake v1)
cal_maxsmd_internal <- function(dat, B = 1200) {
  Tf <- function(zz) max(abs(smd_vec(dat$X, dat$set, zz, dat$sd_pool)))
  Tobs <- Tf(dat$z)
  list(alarm = ref_percentile(Tobs, Tf, function() perm_within(dat$set, dat$z), B), ref = "internal")
}

## 3. matched-set max|SMD| percentile vs whole-pool COMPLETE RANDOMIZATION.
##    EXTERNAL fixed reference. (Jake v2a: complete randomization on the pool)
cal_maxsmd_poolCRE <- function(dat, B = 1500) {
  Tobs <- max(abs(smd_vec(dat$X, dat$set, dat$z, dat$sd_pool)))   # observed: matched-set adjusted
  Tf   <- function(zz) max(abs(smd_vec(dat$X, dat$pool, zz, dat$sd_pool)))  # ref: pool, no strata
  list(alarm = ref_percentile(Tobs, Tf, function() perm_pool(dat$z), B), ref = "external")
}

## 4. fixed-Sigma Mahalanobis of the matched-set mean diff vs whole-pool CRE.
##    EXTERNAL. (Jake v2b: complete randomization, multivariate / covariance metric)
cal_maha_poolCRE <- function(dat, B = 1500) {
  Sig <- cov(dat$X)                                   # FIXED pre-match metric
  Tset <- function(zz) { d <- drop(crossprod(zz - ave(zz, dat$set, FUN = mean), dat$X)); drop(crossprod(d, solve(Sig, d))) }
  Tpool<- function(zz) { d <- drop(crossprod(zz - mean(zz), dat$X)); drop(crossprod(d, solve(Sig, d))) }
  list(alarm = ref_percentile(Tset(dat$z), Tpool, function() perm_pool(dat$z), B), ref = "external")
}

## 5. covariance-ADJUSTED whole-pool CRE: compare the matched design's residual
##    imbalance to a pool CRE whose imbalance is taken AFTER regressing the
##    covariates on the GROUP structure (i.e., the reference also gets credit for
##    the coarse structure adjustment). EXTERNAL. (Jake v2c: CRE + cov adjustment)
cal_maxsmd_poolCRE_adj <- function(dat, B = 1500) {
  ## residualize X on the coarse group means (what adjustment would remove)
  Xr <- apply(dat$X, 2, function(x) x - ave(x, dat$group))
  sdp <- pooled_sd(Xr, dat$z)
  Tobs <- max(abs(smd_vec(Xr, dat$set, dat$z, sdp)))
  Tf   <- function(zz) max(abs(smd_vec(Xr, dat$pool, zz, sdp)))
  list(alarm = ref_percentile(Tobs, Tf, function() perm_pool(dat$z), B), ref = "external")
}

## 6. bare magnitude: observed max|SMD| (fixed SD). NOT a calibration (ref="none").
cal_maxsmd_bare <- function(dat) list(alarm = max(abs(smd_vec(dat$X, dat$set, dat$z, dat$sd_pool))), ref = "none")

## ---- harness: collapse curve (gap FIXED, match tightening) -----------------
collapse_curve <- function(calib, shrinks = c(0.6,0.3,0.15,0.07,0.03,0.012),
                           g = 0.2, reps = 4, seed0 = 1000) {
  rows <- list(); r <- 0
  for (sh in shrinks) {
    vals <- numeric(reps)
    for (k in seq_len(reps)) { dat <- make_matched(sh, g = g, seed = seed0 + 100*match(sh,shrinks) + k)
                               vals[k] <- tryCatch(calib(dat)$alarm, error = function(e) NA_real_) }
    r <- r+1; rows[[r]] <- data.frame(shrink = sh, alarm = mean(vals, na.rm = TRUE))
  }
  out <- do.call(rbind, rows)
  ## collapsed: alarm at the tightest match is near-maximal AND rose from loose->tight
  a_loose <- out$alarm[1]; a_tight <- out$alarm[nrow(out)]
  out_attr <- list(curve = out, ref = tryCatch(calib(make_matched(0.3, g=g, seed=7))$ref, error=function(e) NA),
                   a_loose = a_loose, a_tight = a_tight,
                   collapsed = isTRUE(a_tight >= 0.95 && a_tight - a_loose > 0.1))
  out_attr
}

## ---- harness: magnitude curve (match FIXED, imbalance growing) -------------
## a real calibration's alarm should RISE with the gap; flat alarm = floor-blind.
magnitude_curve <- function(calib, gaps = c(0, 0.1, 0.3, 0.6, 1.0), shrink = 0.1,
                            reps = 4, seed0 = 5000) {
  rows <- list(); r <- 0
  for (g in gaps) {
    vals <- numeric(reps)
    for (k in seq_len(reps)) { dat <- make_matched(shrink, g = g, seed = seed0 + 100*match(g,gaps) + k)
                               vals[k] <- tryCatch(calib(dat)$alarm, error = function(e) NA_real_) }
    r <- r+1; rows[[r]] <- data.frame(gap = g, alarm = mean(vals, na.rm = TRUE))
  }
  do.call(rbind, rows)
}
