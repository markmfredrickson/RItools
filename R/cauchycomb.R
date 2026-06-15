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
##' The statistic is \eqn{T = \sum_k w_k \tan\{(0.5 - p_k)\pi\}} and the combined
##' p-value is \eqn{0.5 - \arctan(T)/\pi}, the survival function of a standard
##' Cauchy.  Validity rests on each \eqn{p_k} being uniform under the null; the
##' Cauchy's stability under averaging is what keeps the combination valid even
##' when the covariates --- and hence their p-values --- are dependent.  A single
##' very small \eqn{p_k} dominates the sum, so the combination is driven by the
##' most imbalanced covariate; large p-values near 1 are guarded so they do not
##' send \code{tan()} to \code{-Inf}.  Missing p-values are dropped.
##'
##' @param p_values numeric vector of p-values to combine.  \code{NA}s are
##'   dropped.
##' @param weights optional non-negative weights, recycled/normalized to sum to
##'   one.  Defaults to equal weights.
##' @return A single combined p-value, or \code{NA_real_} if no non-missing
##'   p-values are supplied.
##' @references Liu, Y. and Xie, J. (2020). Cauchy Combination Test: A Powerful
##'   Test With Analytic p-Value Calculation Under Arbitrary Dependency
##'   Structures. \emph{Journal of the American Statistical Association}
##'   115(529), 393--402.
##' @keywords internal
##' @noRd
acat_pvalue <- function(p_values, weights = NULL) {
  p_values <- p_values[!is.na(p_values)]
  d <- length(p_values)
  if (d == 0L) {
    return(NA_real_)
  }
  if (is.null(weights)) {
    weights <- rep(1 / d, d)
  } else {
    weights <- weights / sum(weights)
  }

  ## Guard the boundaries: tan() diverges to -Inf at p = 1 and +Inf at p = 0.
  ## For tiny p we use the tail approximation tan((0.5 - p) * pi) ~ 1 / (p * pi),
  ## which avoids the precision loss of evaluating tan() right at pi/2 (this is
  ## the same device the published ACAT implementation uses).
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
