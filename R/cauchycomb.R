##' Cauchy combination (ACAT) of per-covariate p-values
##'
##' Combines a set of p-values into a single omnibus p-value using the Cauchy
##' combination test of Liu and Xie (2020).  Unlike the HB08 \eqn{d^2} omnibus,
##' which inverts the (rank-capped) permutation covariance of the adjusted
##' differences, this combination never inverts a covariance.  It therefore has
##' no rank ceiling and stays informative when the number of covariates reaches
##' or exceeds the within-block residual degrees of freedom --- the regime in
##' which \eqn{d^2} freezes at its rank and its p-value parks near 0.46.
##'
##' With \eqn{d} p-values, Liu and Xie's statistic is
##' \eqn{T = d^{-1} \sum_k \tan\{(0.5 - p_k)\pi\}} and the combined
##' p-value is \eqn{0.5 - \arctan(T)/\pi}, the survival function of a standard
##' Cauchy.  Validity rests on each \eqn{p_k} being uniform under the null; the
##' Cauchy's stability under averaging is what keeps the combination valid even
##' when the covariates --- and hence their p-values --- are dependent.  A single
##' very small \eqn{p_k} dominates the sum, so the combination is driven by the
##' most imbalanced covariate.
##'
##' That conversion sends \eqn{p_k = 1} to minus infinity, and in
##' \code{balanceTest()} a covariate's p-value is exactly 1 whenever its adjusted
##' difference is exactly 0, as for a binary covariate with the same treated and
##' control share in every block.  One such covariate would set the combined
##' p-value to 1 whatever the others say.  By default we therefore use the
##' truncated conversion of Gui, Jiang, and Wang (2025), as the riposte package
##' does: with truncation share \eqn{t} each \eqn{p_k} becomes
##' \eqn{\tan\{(0.5 - t p_k)\pi\}}, so \eqn{p_k = 1} becomes
##' \eqn{\tan(-0.4\pi) = -3.08} when \eqn{t = 0.9}, and with \eqn{n} converted
##' values summing to \eqn{S} the combined p-value is
##' \eqn{\min(1, n P(C > S) / t)} for a standard Cauchy \eqn{C}.  At the 0.05
##' level this rejects a true null somewhat more often than Liu and Xie's
##' version (see vignettes/truncated-cauchy-memo.md); at 0.01 it stays close.
##' Missing p-values are dropped.
##'
##' @param p_values numeric vector of p-values to combine.  \code{NA}s are
##'   dropped.
##' @param truncation the share \eqn{t} of the Cauchy distribution kept, in
##'   (0, 1]; 0.9 by default.  \code{truncation = 1} gives Liu and Xie's
##'   untruncated combination.
##' @return A single combined p-value, or \code{NA_real_} if no non-missing
##'   p-values are supplied.
##' @references Liu, Y. and Xie, J. (2020). Cauchy Combination Test: A Powerful
##'   Test With Analytic p-Value Calculation Under Arbitrary Dependency
##'   Structures. \emph{Journal of the American Statistical Association}
##'   115(529), 393--402.
##'
##'   Gui, L., Jiang, Y., and Wang, J. (2025). Aggregating dependent signals
##'   with heavy-tailed combination tests. \emph{Biometrika} 112(4), asaf038.
##' @keywords internal
##' @noRd
acat_pvalue <- function(p_values, truncation = 0.9) {
  if (!is.numeric(truncation) || length(truncation) != 1L || is.na(truncation) ||
      truncation <= 0 || truncation > 1) {
    stop("truncation must be a single number in (0, 1].")
  }
  p_values <- p_values[!is.na(p_values)]
  d <- length(p_values)
  if (d == 0L) {
    return(NA_real_)
  }
  if (truncation < 1) {
    return(truncated_cauchy_log(log(p_values), truncation))
  }

  ## Liu and Xie's untruncated combination, with equal weights.
  ## Guard the boundaries: tan() diverges to -Inf at p = 1 and +Inf at p = 0.
  ## For tiny p we use the tail approximation tan((0.5 - p) * pi) ~ 1 / (p * pi),
  ## which avoids the precision loss of evaluating tan() right at pi/2 (this is
  ## the same device the published ACAT implementation uses).
  weights <- rep(1 / d, d)
  p_values <- pmin(p_values, 1 - .Machine$double.eps)
  is_small <- p_values < 1e-16
  if (any(is_small)) {
    cauchy_stat <- sum((weights[is_small] / p_values[is_small]) / pi) +
      sum(weights[!is_small] * tan((0.5 - p_values[!is_small]) * pi))
  } else {
    cauchy_stat <- sum(weights * tan((0.5 - p_values) * pi))
  }

  0.5 - atan(cauchy_stat) / pi
}

## The truncated combination computed from the logarithms of the p-values, so
## that a p-value too small to convert directly (its converted value overflows)
## still gives a representable combined p-value.  Ported from
## riposte_truncated_cauchy_log() in the riposte package (riposte 0.0.0.9008).
truncated_cauchy_log <- function(log_p, truncation) {
  n <- length(log_p)
  if (any(log_p == -Inf)) return(0)
  ## u = t * p is at most t < 1, so the conversion cot(pi * u) never reaches the
  ## pole at u = 1; below u = 1e-8, cot(pi * u) = 1 / (pi * u) to double precision
  log_u <- log(truncation) + log_p
  small <- log_u < log(1e-8)
  regular <- 1 / tan(pi * exp(log_u[!small]))
  if (!any(small)) {
    S <- sum(regular)
    return(min(1, n * stats::pcauchy(S, lower.tail = FALSE) / truncation))
  }
  ## at least one converted value exceeds 3e7: add the large ones on the log
  ## scale, and use the Cauchy tail 1 / (pi * S), exact to double precision there
  log_big <- -log_u[small] - log(pi)
  top <- max(log_big)
  log_S <- top + log(sum(exp(log_big - top)) + sum(regular) * exp(-top))
  log_p_comb <- log(n) - log_S - log(pi) - log(truncation)
  min(1, exp(log_p_comb))
}
