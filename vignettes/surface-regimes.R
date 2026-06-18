## ---------------------------------------------------------------------------
## surface-regimes.R
##
## Characterize Jake's 2D surface (max|SMD|, omnibus p) traced by RE-ORGANIZING
## fixed units, across regimes and ladder types, to answer three questions:
##
##  Q1 (race)   : the omnibus chisq = M x P, M = multivariate SIZE
##                (dbar' Sigma_pool^{-1} dbar), P = PRECISION (Sigma_pool / V_d
##                summarized as chisq / M). Reorganizing trades M against P.
##                Whether refinement makes p RISE or FALL is whoever wins.
##  Q2 (regimes): if the organizing score is a GOOD proxy for the confounder,
##                M falls fast and p rises (balance wins); if a WEAK proxy, M
##                floors while P climbs and p falls (the collapse-like illusion).
##  Q3 (jitter) : is the non-monotonicity deep, or an artifact of comparing
##                NON-NESTED partitions? Compare a non-nested quantile-cut ladder
##                with a genuinely NESTED (agglomerative) ladder.
##
## Output: printed summaries + RDS artifacts under vignettes/ for the
## interrogation workflow to read and adversarially verify.
## ---------------------------------------------------------------------------

.pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."
suppressMessages(devtools::load_all(.pkg_root, export_all = FALSE, quiet = TRUE))

## ---- fixed-unit DGP -------------------------------------------------------
## N per arm, K correlated covariates (equicorrelation rho). Treated shifted by
## `delta` along a fixed unit direction u. The organizing score is a proxy for
## the position along u, with controllable noise: small noise = STRONG proxy
## (refinement removes the confounder), large noise = WEAK proxy (residual bias
## persists). delta = 0 is the null.
make_units <- function(seed, N = 200, K = 4, rho = 0.5, delta = 0.5,
                       score_noise = 0.6) {
  set.seed(seed)
  R <- matrix(rho, K, K); diag(R) <- 1
  L <- t(chol(R))
  draw <- function(n) t(L %*% matrix(rnorm(K * n), nrow = K))
  u  <- rep(1, K) / sqrt(K)
  Xc <- draw(N)
  Xt <- draw(N) + matrix(delta * u, N, K, byrow = TRUE)
  X  <- rbind(Xt, Xc); colnames(X) <- paste0("x", seq_len(K))
  z  <- c(rep(1L, N), rep(0L, N))
  s  <- as.numeric(X %*% u) + rnorm(2 * N, sd = score_noise)
  list(df = data.frame(z = z, score = s, X), covs = colnames(X))
}

## ---- two ladders ----------------------------------------------------------
## non-nested: quantile cut of the score into S bins (boundaries move with S)
ladder_cut <- function(df, S) {
  if (S <= 1) return(rep(1L, nrow(df)))
  br <- quantile(df$score, seq(0, 1, length.out = S + 1))
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(df$score, breaks = unique(br), labels = FALSE))
}
## nested: agglomerative clustering of the score; cutree(k) refines cutree(k-1)
make_hc <- function(df) hclust(dist(df$score), method = "ward.D2")
ladder_nested <- function(hc, S) if (S <= 1) rep(1L, length(hc$order)) else cutree(hc, k = S)

## ---- evaluate one organization -------------------------------------------
eval_org <- function(df, covs, fac, Sigma_pool) {
  d2 <- df; d2$.strat <- factor(fac)
  tab <- table(d2$.strat, d2$z)
  ok  <- rownames(tab)[tab[, 1] > 0 & tab[, 2] > 0]
  d2  <- d2[as.character(d2$.strat) %in% ok, , drop = FALSE]
  d2$.strat <- droplevels(d2$.strat)
  if (nlevels(d2$.strat) < 1 || nrow(d2) < 4) return(NULL)
  tryCatch({
    bt  <- balanceTest(reformulate(c(covs, "strata(.strat)"), response = "z"), data = d2)
    smd <- bt$results[, "std.diff", ".strat"]
    adj <- bt$results[, "adj.diff", ".strat"]
    ov  <- bt$overall[".strat", ]
    M_maha <- drop(crossprod(adj, solve(Sigma_pool, adj)))
    tb <- table(d2$.strat, d2$z)
    ess <- sum(2 * tb[, 1] * tb[, 2] / (tb[, 1] + tb[, 2]))
    data.frame(n_strata = nlevels(d2$.strat), n_units = nrow(d2), ess = ess,
               max_smd = max(abs(smd)), mean_smd = mean(abs(smd)),
               M_diag = sum(smd^2), M_maha = M_maha,
               chisq = ov[["chisquare"]], df = ov[["df"]], p = ov[["p.value"]],
               P = chisq_over_M(ov[["chisquare"]], M_maha))
  }, error = function(e) NULL)
}
chisq_over_M <- function(chisq, M) if (is.na(M) || M <= 0) NA_real_ else chisq / M

## ---- sweep one dataset over one ladder ------------------------------------
sweep_one <- function(seed, regime, ladder, S_grid) {
  u  <- make_units(seed, delta = regime$delta, score_noise = regime$score_noise)
  df <- u$df; covs <- u$covs
  Sigma_pool <- cov(as.matrix(df[, covs]))
  hc <- if (ladder == "nested") make_hc(df) else NULL
  rows <- lapply(S_grid, function(S) {
    fac <- if (ladder == "nested") ladder_nested(hc, S) else ladder_cut(df, S)
    r <- eval_org(df, covs, fac, Sigma_pool)
    if (!is.null(r)) { r$S_req <- S; r$seed <- seed; r$ladder <- ladder; r$regime <- regime$name }
    r
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

## ---- experiment grid ------------------------------------------------------
regimes <- list(
  strong = list(name = "strong", delta = 0.5, score_noise = 0.15),  # good proxy
  weak   = list(name = "weak",   delta = 0.5, score_noise = 1.50),  # weak proxy
  null   = list(name = "null",   delta = 0.0, score_noise = 0.60)   # no confounder
)
S_grid  <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 26, 34, 44, 56, 72, 92, 120, 150)
seeds   <- 1:30

all <- list()
for (rg in regimes) for (ld in c("cut", "nested")) {
  key <- paste(rg$name, ld, sep = "_")
  cat("running", key, "...\n")
  all[[key]] <- do.call(rbind, lapply(seeds, function(s) sweep_one(s, rg, ld, S_grid)))
}
dat <- do.call(rbind, all); rownames(dat) <- NULL
saveRDS(dat, "vignettes/surface-regimes.rds")

## ---- seed-averaged trend (separate structure from dataset jitter) ---------
agg <- aggregate(cbind(max_smd, M_diag, M_maha, chisq, P, p) ~ regime + ladder + S_req,
                 data = dat, FUN = median)
agg <- agg[order(agg$regime, agg$ladder, agg$S_req), ]

## within-dataset jitter: average over seeds of the #direction-reversals in p
## along the ladder (a clean structural-jitter measure, NOT seed noise since p
## is the analytic chisq p-value for a fixed dataset)
reversals <- function(p) { d <- sign(diff(p)); d <- d[d != 0]; if (length(d) < 2) 0 else sum(diff(d) != 0) }
jit <- aggregate(p ~ regime + ladder + seed, data = dat[order(dat$S_req), ],
                 FUN = reversals)
jit_summary <- aggregate(p ~ regime + ladder, data = jit, FUN = mean)
names(jit_summary)[3] <- "mean_p_reversals"

cat("\n================ seed-median trend (by regime x ladder) ================\n")
for (rg in names(regimes)) for (ld in c("cut", "nested")) {
  sub <- agg[agg$regime == rg & agg$ladder == ld, ]
  cat(sprintf("\n--- regime=%s  ladder=%s ---\n", rg, ld))
  pr <- sub[, c("S_req", "max_smd", "M_maha", "P", "chisq", "p")]
  pr$max_smd <- round(pr$max_smd, 3); pr$M_maha <- round(pr$M_maha, 4)
  pr$P <- round(pr$P); pr$chisq <- round(pr$chisq, 1); pr$p <- signif(pr$p, 3)
  print(pr, row.names = FALSE)
}

cat("\n================ structural jitter: mean # p-reversals along ladder ====\n")
print(jit_summary, row.names = FALSE, digits = 3)

## ---- direction-of-refinement summary --------------------------------------
## does p go UP or DOWN from coarsest to finest (seed-median)?
cat("\n================ refinement direction (coarse -> fine, seed-median) ====\n")
for (rg in names(regimes)) for (ld in c("cut", "nested")) {
  sub <- agg[agg$regime == rg & agg$ladder == ld, ]
  sub <- sub[order(sub$S_req), ]
  cat(sprintf("regime=%-6s ladder=%-6s  max_smd %.3f->%.3f   p %.3g->%.3g   M_maha %.3f->%.3f   P %.0f->%.0f\n",
              rg, ld, sub$max_smd[1], sub$max_smd[nrow(sub)],
              sub$p[1], sub$p[nrow(sub)], sub$M_maha[1], sub$M_maha[nrow(sub)],
              sub$P[1], sub$P[nrow(sub)]))
}

## ---- plots ----------------------------------------------------------------
png("vignettes/surface-regimes.png", width = 2000, height = 1500, res = 175)
op <- par(mfrow = c(2, 3), mar = c(4.2, 4.4, 2.8, 1))
cols <- c(strong = "#1b9e77", weak = "#d95f02", null = "#7570b3")
for (ld in c("cut", "nested")) {
  ## p vs strata
  plot(NA, xlim = range(S_grid), ylim = c(1e-12, 1), log = "xy",
       xlab = "number of strata (finer ->)", ylab = "omnibus p (d^2)",
       main = sprintf("p vs refinement [%s ladder]", ld))
  abline(h = 0.05, lty = 3, col = "grey60")
  for (rg in names(regimes)) {
    sub <- agg[agg$regime == rg & agg$ladder == ld, ]; sub <- sub[order(sub$S_req), ]
    lines(sub$S_req, pmax(sub$p, 1e-12), type = "b", pch = 19, col = cols[rg])
  }
  legend("bottomleft", names(regimes), col = cols, lwd = 2, pch = 19, bty = "n")
  ## max_smd vs strata
  plot(NA, xlim = range(S_grid), ylim = c(0, max(agg$max_smd)), log = "x",
       xlab = "number of strata (finer ->)", ylab = "max |SMD|",
       main = sprintf("max|SMD| vs refinement [%s]", ld))
  for (rg in names(regimes)) {
    sub <- agg[agg$regime == rg & agg$ladder == ld, ]; sub <- sub[order(sub$S_req), ]
    lines(sub$S_req, sub$max_smd, type = "b", pch = 19, col = cols[rg])
  }
  ## the surface itself
  plot(NA, xlim = c(0, max(agg$max_smd)), ylim = c(1e-12, 1), log = "y",
       xlab = "max |SMD| (size)", ylab = "omnibus p (precision-laden)",
       main = sprintf("the surface [%s]", ld))
  abline(h = 0.05, lty = 3, col = "grey60")
  for (rg in names(regimes)) {
    sub <- agg[agg$regime == rg & agg$ladder == ld, ]; sub <- sub[order(sub$S_req), ]
    points(sub$max_smd, pmax(sub$p, 1e-12), type = "b", pch = 19, col = cols[rg])
  }
}
par(op); dev.off()
cat("\nplot written to vignettes/surface-regimes.png\n")
cat("data written to vignettes/surface-regimes.rds\n")
