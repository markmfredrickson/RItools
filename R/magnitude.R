## Magnitude balance report -- two numbers that, unlike the omnibus p-value, do
## not collapse as the match tightens.  Form (a) is a descriptive reader on the
## standardized differences balanceTest() already computes; form (b) is an
## EXPERIMENTAL, internal structure-adjusted whole-pool randomization percentile.
## Neither touches the omnibus: the magnitude report is reported ALONGSIDE the
## d^2 test, never as a pre-filter on the covariates fed to it (decision D).

##' Descriptive magnitude balance report
##'
##' Reads the per-covariate standardized differences that
##' \code{\link{balanceTest}} already computes (\code{std.diff = adj.diff /
##' pooled.sd}) and presents them as a magnitude report: each covariate's
##' standardized difference and its absolute value, whether it exceeds a settable
##' threshold, and the global \code{max|SMD|}.  This is a fixed-scale magnitude --
##' it does not collapse as the match tightens, unlike the omnibus p-value -- and
##' it is meant to be read alongside the omnibus, not used to pre-filter the
##' covariates fed to it.
##'
##' The default \code{threshold} of 0.25 is a matching-literature convention (see
##' \code{vignettes/balance-threshold-provenance.md}); it is not a law of nature.
##' A substantively chosen, design-stage tolerance -- for example, a maximum age
##' difference of so many years -- is usually a better yardstick than a universal
##' standardized cutoff.  Note also that the denominator here is the pooled
##' within-group SD that RItools uses; some traditions standardize by the
##' treated-group SD.
##'
##' @param object a \code{balanceTest()} (or \code{xBalance()}) result; only its
##'   \code{$results} table is read.
##' @param threshold numeric; absolute standardized differences above this are
##'   marked as exceedances.  Default 0.25.
##' @param stratification character vector of stratification names to report
##'   (a subset of \code{dimnames(object$results)[["strata"]]}, which includes
##'   the unstratified \code{"--"} comparison).  \code{NULL} (default) reports all.
##' @return an object of class \code{balanceMagnitude}: a named list of
##'   data frames, one per stratification, each with columns \code{variable},
##'   \code{std.diff}, \code{abs.std.diff}, and \code{exceeds}, sorted by
##'   \code{abs.std.diff}.  Each data frame carries \code{max.abs.std.diff} and
##'   \code{threshold} attributes.
##' @seealso \code{\link{balanceTest}}
##' @examples
##' data(nuclearplants, package = "RItools")
##' bt <- balanceTest(pr ~ date + t1 + t2 + cap, data = nuclearplants)
##' balanceMagnitude(bt, threshold = 0.25)
##' @export
balanceMagnitude <- function(object, threshold = 0.25, stratification = NULL) {
  res <- object[["results"]]
  if (is.null(res))
    stop("object has no $results table; pass a balanceTest()/xBalance() result")
  if (!"std.diff" %in% dimnames(res)[["stat"]])
    stop("object$results has no std.diff column")
  strata_all <- dimnames(res)[["strata"]]
  if (is.null(stratification)) {
    stratification <- strata_all
  } else {
    bad <- setdiff(stratification, strata_all)
    if (length(bad))
      stop("unknown stratification(s): ", paste(bad, collapse = ", "),
           "; available: ", paste(strata_all, collapse = ", "))
  }
  varnames <- dimnames(res)[["vars"]]
  tables <- lapply(stratification, function(st) {
    sd <- as.numeric(res[, "std.diff", st])
    a  <- abs(sd)
    df <- data.frame(variable = varnames, std.diff = sd, abs.std.diff = a,
                     exceeds = a > threshold, stringsAsFactors = FALSE)
    df <- df[order(-df$abs.std.diff), , drop = FALSE]
    rownames(df) <- NULL
    attr(df, "max.abs.std.diff") <- if (all(is.na(a))) NA_real_ else max(a, na.rm = TRUE)
    attr(df, "threshold")        <- threshold
    attr(df, "stratification")   <- st
    df
  })
  names(tables) <- stratification
  structure(tables, class = "balanceMagnitude", threshold = threshold)
}

##' Print a magnitude balance report
##'
##' @param x a \code{balanceMagnitude} object
##' @param digits number of digits for the standardized differences
##' @param ... ignored
##' @return \code{x}, invisibly
##' @export
print.balanceMagnitude <- function(x, digits = 3, ...) {
  for (st in names(x)) {
    df <- x[[st]]
    cat(sprintf("Stratification: %s\n", st))
    mark <- ifelse(is.na(df$exceeds), "?", ifelse(df$exceeds, "*", ""))
    out <- data.frame(variable = df$variable,
                      std.diff = round(df$std.diff, digits),
                      abs.std.diff = round(df$abs.std.diff, digits),
                      " " = mark, check.names = FALSE, stringsAsFactors = FALSE)
    print(out, row.names = FALSE)
    cat(sprintf("max|SMD| = %s  (threshold %s; * exceeds)\n\n",
                format(round(attr(df, "max.abs.std.diff"), digits)),
                format(attr(df, "threshold"))))
  }
  invisible(x)
}

## ---------------------------------------------------------------------------
## Form (b): structure-adjusted whole-pool complete-randomization percentile.
## EXPERIMENTAL and INTERNAL.  Productionizes the finding in
## vignettes/impossibility-pressure-test-memo.md (sec 4): a fixed external
## reference (complete randomization on the whole pool) escapes the collapse but
## is floor-blind unless the covariates are first adjusted for the coarse
## confounding structure.  The adjustment is recovered UNSUPERVISED from the
## observed covariates by k-means; residualizing each covariate on the OTHERS
## (the naive move) is deliberately avoided, because it absorbs the
## common-direction gap and kills magnitude sensitivity.
## ---------------------------------------------------------------------------

## residualize each covariate on a k-means clustering of the (scaled) covariates
.residualize_on_kmeans <- function(X, k) {
  X <- as.matrix(X)
  km <- tryCatch(stats::kmeans(scale(X), centers = k, nstart = 10L, iter.max = 50L),
                 error = function(e) NULL)
  if (is.null(km)) return(X)
  glab <- factor(km$cluster)
  apply(X, 2, function(col) col - ave(col, glab))
}

## within-cluster dispersion W_k for the gap statistic: the pooled within-cluster
## sum of squares (kmeans tot.withinss); for k = 1 it is the total sum of squares.
.cluster_dispersion <- function(Xs, k) {
  if (k <= 1L) return(sum(scale(Xs, scale = FALSE)^2))
  km <- tryCatch(stats::kmeans(Xs, centers = k, nstart = 5L, iter.max = 50L),
                 error = function(e) NULL)
  if (is.null(km)) return(NA_real_)
  km$tot.withinss
}

## Data-driven number of coarse groups by the gap statistic (Tibshirani, Walther,
## Hastie 2001).  We compare log W_k for the data to its expectation under a
## uniform reference over the bounding box of the standardized covariates, and
## take the smallest k whose gap is within one standard error of the next k's.
## This recovers the coarse confounding structure the matching exploited without
## us having to assert how many groups there are, which is why k is chosen here
## rather than fixed.  B_gap is kept small because this sits inside an already
## Monte-Carlo percentile.
.choose_k_gap <- function(X, kmax = NULL, B_gap = 15L) {
  Xs <- scale(as.matrix(X))
  n <- nrow(Xs); p <- ncol(Xs)
  if (is.null(kmax)) kmax <- max(2L, min(8L, floor(n / 10)))
  ks <- 1:kmax
  logWk <- vapply(ks, function(k) log(.cluster_dispersion(Xs, k)), numeric(1))
  rng <- apply(Xs, 2, range)
  logWref <- matrix(NA_real_, B_gap, length(ks))
  for (b in seq_len(B_gap)) {
    Xb <- vapply(seq_len(p),
                 function(j) stats::runif(n, rng[1, j], rng[2, j]), numeric(n))
    logWref[b, ] <- vapply(ks, function(k) log(.cluster_dispersion(Xb, k)), numeric(1))
  }
  gap <- colMeans(logWref) - logWk
  s   <- apply(logWref, 2, stats::sd) * sqrt(1 + 1 / B_gap)
  for (i in seq_len(length(ks) - 1L)) {
    if (is.finite(gap[i]) && is.finite(gap[i + 1L]) &&
        gap[i] >= gap[i + 1L] - s[i + 1L]) return(ks[i])
  }
  ks[length(ks)]
}

## pooled within-group SD per covariate (the magnitude denominator)
.magnitude_pooled_sd <- function(X, z) {
  X <- as.matrix(X)
  apply(X, 2, function(x) {
    v1 <- var(x[z == 1]); v0 <- var(x[z == 0])
    n1 <- sum(z == 1); n0 <- sum(z == 0)
    sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
  })
}

## weighted within-stratum standardized mean difference, max over covariates
.within_stratum_max_smd <- function(X, strata, z, sd_pool) {
  X <- as.matrix(X); K <- ncol(X); num <- numeric(K); wsum <- 0
  for (idx in split(seq_len(nrow(X)), strata)) {
    zz <- z[idx]; n1 <- sum(zz); n0 <- length(zz) - n1
    if (n1 == 0 || n0 == 0) next
    w  <- n1 * n0 / (n1 + n0)
    mt <- colMeans(X[idx, , drop = FALSE][zz == 1, , drop = FALSE])
    mc <- colMeans(X[idx, , drop = FALSE][zz == 0, , drop = FALSE])
    num <- num + w * (mt - mc); wsum <- wsum + w
  }
  max(abs((num / wsum) / sd_pool))
}

##' Structure-adjusted whole-pool randomization percentile (EXPERIMENTAL)
##'
##' EXPERIMENTAL and not exported.  Computes a calibrated balance number that,
##' unlike the omnibus p-value, does not collapse as the match tightens: it
##' compares the matched design's residual \code{max|SMD|} to the distribution of
##' \code{max|SMD|} under complete randomization on the whole pool, after
##' residualizing the covariates on a coarse structure recovered by k-means.  See
##' \code{vignettes/impossibility-pressure-test-memo.md}.
##'
##' The number of recovered groups \code{k} is chosen data-driven by default (the
##' gap statistic of Tibshirani, Walther and Hastie, 2001); pass an integer to fix
##' it.  \code{centering} selects the contender.  \code{"conservative"} (default)
##' compares a matched-set-centered observed statistic to a pool-centered
##' reference: a one-sided "better than a whole-pool coin flip" check that is
##' downward biased, not an exactly uniform tail probability.  \code{"uniform"}
##' centers both observed and reference on the pool, giving an exactly uniform
##' percentile under complete randomization, at the cost of discarding the
##' matched-set structure in the observed statistic.  The denominator is the
##' pooled within-group SD.  Defaults that may still want Jake/Ben review: the
##' gap-statistic range and reference-draw count, and the number of
##' complete-randomization draws \code{B}.  See
##' \code{vignettes/impossibility-pressure-test-memo.md}.
##'
##' @param X numeric matrix of covariates (units in rows)
##' @param z 0/1 treatment indicator
##' @param strata factor of matched-set membership
##' @param k number of coarse groups to recover by k-means; \code{NULL} (default)
##'   chooses it data-driven via the gap statistic
##' @param centering \code{"conservative"} (matched-set observed vs pool reference;
##'   default) or \code{"uniform"} (pool vs pool); see Details
##' @param B number of complete-randomization draws
##' @param B_gap number of uniform reference datasets for the gap statistic
##' @return upper-tail percentile in \code{[0,1]}
##' @keywords internal
poolCRE_adjusted_percentile <- function(X, z, strata, k = NULL,
                                        centering = c("conservative", "uniform"),
                                        B = 1000L, B_gap = 15L) {
  centering <- match.arg(centering)
  X <- as.matrix(X); z <- as.integer(z); strata <- factor(strata)
  if (is.null(k)) k <- .choose_k_gap(X, B_gap = B_gap)
  Xr   <- .residualize_on_kmeans(X, k)
  sdp  <- .magnitude_pooled_sd(Xr, z)
  pool <- factor(rep(1L, length(z)))
  ## conservative keeps the matched-set structure in the observed statistic
  ## (smaller, downward biased); uniform centers the observed on the pool like the
  ## reference, so the percentile is exactly uniform under complete randomization.
  obs_strata <- if (centering == "conservative") strata else pool
  Tobs  <- .within_stratum_max_smd(Xr, obs_strata, z, sdp)
  Tnull <- replicate(B, .within_stratum_max_smd(Xr, pool, sample(z), sdp))
  mean(Tnull <= Tobs - 1e-12)
}
