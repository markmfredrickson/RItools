################################################################################
## Sigma_x balance test:  T_new = d' Sigma_x^{-1} d
##
## Internal functions for an alternative omnibus statistic that standardizes
## the adjusted differences vector by an external covariance Sigma_x rather
## than by the permutation covariance Cov(d) used by the Hansen-Bowers d^2
## statistic in HB08.  See ?balanceTest, argument `sigma_x_test`.
##
## Layered API (all unexported):
##   sigma_x_pvalue(q, lambda, EM, VarM, method)        scalar p-value
##   default_sigma_x(X, strata)                          within-stratum-pooled cov
##   randomization_cov_d(X, strata, n1_per_stratum)      closed-form Cov(d_raw)
##   sigma_x_T_moments(X, strata, n1_per_stratum, sigma_x)
##                                                       exact E[T] and Var[T]
##                                                       under within-stratum SRSWOR
##   sigma_x_test(X, strata, z, sigma_x, null, n_simulate)
##                                                       main entry on raw inputs
##   sigma_x_inferentials(aggDesign, stratification_name, sigma_x, null,
##                        n_simulate)
##                                                       hook used by balanceTest()
##
## The moment machinery in sigma_x_T_moments() is a port of work from the
## i113-highermoments branch (Mark Fredrickson) -- see comments where the
## Finucan-style coefficient formulas are reused.
################################################################################

##' P-value for a quadratic form in finite-dimensional Normals
##'
##' Computes \eqn{P(T > q)} for \eqn{T = \sum_k \lambda_k \chi^2_1} via one of:
##' \itemize{
##'   \item \code{"satterthwaite_asymptotic"}: a Satterthwaite \eqn{a\chi^2_\nu}
##'         moment match using the Gaussian identities \eqn{E[T] = \sum\lambda_k}
##'         and \eqn{Var[T] = 2\sum\lambda_k^2}.  Asymptotic in \eqn{n}.
##'   \item \code{"satterthwaite_finite"}: the same moment match, but with the
##'         caller supplying \eqn{E[T]} and \eqn{Var[T]} computed under the
##'         exact randomization distribution.  Honest at finite \eqn{n}.
##'   \item \code{"imhof"}, \code{"davies"}: analytic CDF inversion via
##'         \pkg{CompQuadForm}.  Also assume Gaussian \eqn{d}.
##' }
##'
##' @param q observed test statistic value
##' @param lambda eigenvalues of \eqn{\Sigma_x^{-1/2} V_d \Sigma_x^{-1/2}};
##'   not needed for \code{satterthwaite_finite}
##' @param EM scalar \eqn{E[T]}; required for \code{satterthwaite_finite}
##' @param VarM scalar \eqn{Var[T]}; required for \code{satterthwaite_finite}
##' @param method one of the four backends above
##' @keywords internal
sigma_x_pvalue <- function(q,
                            lambda = NULL,
                            EM = NULL,
                            VarM = NULL,
                            method = c("satterthwaite_asymptotic",
                                       "satterthwaite_finite",
                                       "imhof", "davies")) {
  method <- match.arg(method)

  if (method == "satterthwaite_finite") {
    if (is.null(EM) || is.null(VarM)) {
      stop("'satterthwaite_finite' requires both EM and VarM")
    }
    return(satterthwaite_pvalue(q, EM, VarM))
  }

  ## Asymptotic backends all need the eigenvalue spectrum
  if (is.null(lambda)) {
    stop("method '", method, "' requires lambda")
  }
  ## Drop numerical zeros so the Satterthwaite/Imhof formulas don't see
  ## meaningless tiny eigenvalues
  if (length(lambda) > 0) {
    keep <- lambda > sqrt(.Machine$double.eps) * max(lambda, 1)
    lambda <- lambda[keep]
  }
  if (length(lambda) == 0) return(1)  # degenerate; null cannot reject

  if (method == "satterthwaite_asymptotic") {
    return(satterthwaite_pvalue(q,
                                 EM   = sum(lambda),
                                 VarM = 2 * sum(lambda^2)))
  }
  if (method %in% c("imhof", "davies")) {
    if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
      stop("method '", method, "' requires the CompQuadForm package; ",
           "install it or pick a different null backend")
    }
    fn <- if (method == "imhof") CompQuadForm::imhof else CompQuadForm::davies
    return(fn(q, lambda = lambda)$Qq)
  }
  stop("unreachable")
}

## Helper: Satterthwaite a*chi^2_v p-value from (E[T], Var[T]).
##   E[a chi^2_v]   = a v   = EM
##   Var[a chi^2_v] = 2 a^2 v = VarM
## Solving for a and v:  a = VarM / (2 EM),  v = 2 EM^2 / VarM,
## then  P(T > q) = 1 - F_chisq(q / a, df = v).
satterthwaite_pvalue <- function(q, EM, VarM) {
  if (!is.finite(EM) || !is.finite(VarM) || EM <= 0 || VarM <= 0) {
    return(NA_real_)
  }
  a <- VarM / (2 * EM)
  v <- 2 * EM^2 / VarM
  pchisq(q / a, df = v, lower.tail = FALSE)
}

##' Within-stratum-pooled sample covariance, the proposed default for Sigma_x
##'
##' \deqn{\Sigma_x = \frac{1}{N - K} \sum_{s} \sum_{i \in s}
##'                  (x_i - \bar x_s)(x_i - \bar x_s)'}{(1/(N-K)) sum_s sum_{i in s} (x_i - xbar_s)(x_i - xbar_s)'}
##'
##' This is the multivariate analogue of the pooled-SD denominator the existing
##' \code{balanceTest()} already uses for the per-covariate standardized
##' differences.  In the unstratified case it reduces to \code{cov(X)}.
##' @keywords internal
default_sigma_x <- function(X, strata) {
  X <- as.matrix(X)
  K <- nlevels(strata)
  N <- nrow(X)
  if (N - K < 1) {
    stop("need at least one unit beyond the number of strata")
  }
  by_s <- split(seq_len(N), strata)
  S <- matrix(0, ncol(X), ncol(X))
  for (idx in by_s) {
    if (length(idx) < 2) next
    Xc <- scale(X[idx, , drop = FALSE], scale = FALSE)
    S <- S + crossprod(Xc)
  }
  out <- S / (N - K)
  dimnames(out) <- list(colnames(X), colnames(X))
  out
}

##' Closed form for Cov((z - pi)' X) under within-stratum SRSWOR
##'
##' Under stratified simple random sampling without replacement, the second-
##' order moments of \eqn{z} are
##' \eqn{Var(z_i) = \pi_s(1-\pi_s)} and
##' \eqn{Cov(z_i, z_j) = -\pi_s(1-\pi_s)/(n_s - 1)} for \eqn{i,j} in the same
##' stratum, and zero across strata.  Substituting and using the within-
##' stratum centering of \eqn{(z - \pi)} gives
##' \deqn{Cov(d_{raw}) = \sum_s \frac{n_{1s} n_{0s}}{n_s (n_s - 1)} S_{xs}}
##' where \eqn{S_{xs} = \sum_{i \in s} (x_i - \bar x_s)(x_i - \bar x_s)'}.
##' @keywords internal
randomization_cov_d <- function(X, strata, n1_per_stratum) {
  X <- as.matrix(X)
  K <- nlevels(strata)
  if (length(n1_per_stratum) != K) {
    stop("n1_per_stratum has wrong length")
  }
  n_per_s <- tabulate(as.integer(strata))
  by_s <- split(seq_len(nrow(X)), strata)
  V <- matrix(0, ncol(X), ncol(X))
  for (s in seq_len(K)) {
    idx  <- by_s[[s]]
    n_s  <- n_per_s[s]
    n1_s <- n1_per_stratum[s]
    n0_s <- n_s - n1_s
    if (n_s < 2) next
    Xc <- scale(X[idx, , drop = FALSE], scale = FALSE)
    V  <- V + (n1_s * n0_s / (n_s * (n_s - 1))) * crossprod(Xc)
  }
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

##' Exact \eqn{E[T]} and \eqn{Var[T]} under within-stratum SRSWOR
##'
##' Implements the moment pair the finite-sample Satterthwaite backend needs.
##'
##' Strategy: rotate \eqn{X \to Y = X \Sigma_x^{-1/2}} so that the test
##' statistic in the rotated frame is \eqn{\|d_y\|^2}, then apply per-stratum
##' first- and second-order moment formulas to compute \eqn{E[\|d_y\|^2]}
##' (= \eqn{trace(Cov(d_y))}) and \eqn{Var[\|d_y\|^2]}.  For independent strata,
##' \deqn{Var(\|T\|^2) = \sum_s Var(\|T_s\|^2) +
##'        4 \sum_{i<j} \sum_{a,b} Cov(T_i)_{ab} Cov(T_j)_{ab}}
##' Per-stratum \eqn{Var(\|T_s\|^2)} comes from the Finucan-style fourth-order
##' coefficient formulas ported from the i113-highermoments branch.  At
##' \eqn{n_s \le 3} the formulas degenerate; we follow the branch's special
##' case (collapse to the empirical fourth moment), which gives the correct
##' answer of zero for pair-matched (\eqn{n_s = 2}) strata.
##' @keywords internal
sigma_x_T_moments <- function(X, strata, n1_per_stratum, sigma_x) {
  X <- as.matrix(X)
  ## Rotate by Sigma_x^{-1/2} so the test stat becomes ||d_y||^2 in the
  ## rotated frame.  XtX_pseudoinv_sqrt returns a p x r matrix with
  ## r = rank(sigma_x), so this works even when Sigma_x is rank-deficient.
  rot <- XtX_pseudoinv_sqrt(sigma_x, mat.is.XtX = TRUE)
  Y   <- as.matrix(X %*% rot)
  r   <- ncol(Y)

  ## E[T] = trace(Cov(d_y)) = sum of per-stratum traces
  V_y <- randomization_cov_d(Y, strata, n1_per_stratum)
  EM  <- sum(diag(V_y))

  ## Per-stratum first- and second-order moment arrays for Var(||d_y||^2)
  K <- nlevels(strata)
  by_s <- split(seq_len(nrow(Y)), strata)
  n_per_s <- tabulate(as.integer(strata))

  first_order  <- array(0, c(K, r, r))   # Cov(T_s)
  second_order <- array(0, c(K, r, r))   # Cov(T_{s,a}^2, T_{s,b}^2)

  for (s in seq_len(K)) {
    idx  <- by_s[[s]]
    n_s  <- n_per_s[s]
    n1_s <- n1_per_stratum[s]
    n0_s <- n_s - n1_s
    n1n0 <- n1_s * n0_s
    if (n_s < 2) next

    Yc <- scale(Y[idx, , drop = FALSE], scale = FALSE)  # within-stratum centered

    ## First-order: Cov(T_s) = (n1n0 / (n_s (n_s - 1))) * S_ys
    first_order[s, , ] <- (n1n0 / (n_s * (n_s - 1))) * crossprod(Yc)

    ## Second-order: per-stratum Cov(T_{s,a}^2, T_{s,b}^2) via Finucan.  See
    ## i113-highermoments R/Stratified.R::strata_t2_covariance_matrices().
    mu11    <- crossprod(Yc) / n_s                 # (1/n_s) sum (y - ybar)(y - ybar)'
    mu22    <- crossprod(Yc^2) / n_s               # 4th-order monomial mean
    mu2     <- colSums(Yc^2) / n_s                 # diagonal of mu11
    mu2_mu2 <- tcrossprod(mu2)                     # E[T_a^2] E[T_b^2] for ind. moments

    coef1   <- n1n0 / (n_s - 1)
    if (n_s > 3) {
      coef2  <- coef1 / ((n_s - 2) * (n_s - 3))
      coef2a <- n_s * (n_s + 1) - 6 * n1n0
      coef2b <- n_s * (n_s - 1 - n1n0)
      mean_22 <- coef2 * (coef2a * mu22 - coef2b * (2 * mu11^2 + mu2_mu2))
    } else {
      ## Small-stratum branch (n_s = 2 or 3).  At n_s = 2 (pair match) the
      ## per-stratum second-order array collapses to zero, matching the fact
      ## that ||T_s||^2 is deterministic for pair matches.
      mean_22 <- mu22
    }
    mean2_mean2 <- coef1^2 * mu2_mu2
    second_order[s, , ] <- mean_22 - mean2_mean2
  }

  ## Sum of per-stratum Var(||T_s||^2) = sum over s of sum_{a,b} Cov(T_a^2, T_b^2)
  diag_term <- sum(second_order)

  ## Cross-stratum: 4 * sum_{i<j} sum(Cov(T_i) o Cov(T_j))  (element-wise)
  cross_term <- 0
  if (K >= 2) {
    for (i in seq_len(K - 1)) {
      for (j in seq.int(i + 1, K)) {
        cross_term <- cross_term + sum(first_order[i, , ] * first_order[j, , ])
      }
    }
    cross_term <- 4 * cross_term
  }
  VarM <- diag_term + cross_term

  list(EM = EM, VarM = VarM, rank = r)
}

##' The Sigma_x balance test from raw inputs
##'
##' Mid-level entry point.  Given a raw covariate matrix, a strata factor, a
##' treatment vector, an optional Sigma_x, and a choice of null backend,
##' computes the statistic \eqn{T = d' \Sigma_x^{-1} d} and its p-value.
##' @keywords internal
sigma_x_test <- function(X, strata, z,
                          sigma_x = NULL,
                          null = c("satterthwaite_finite",
                                   "satterthwaite_asymptotic",
                                   "imhof", "davies", "simulate"),
                          n_simulate = 1000) {
  null <- match.arg(null)
  X <- as.matrix(X)
  stopifnot(is.factor(strata),
            length(strata) == nrow(X),
            length(z) == nrow(X))

  if (is.null(sigma_x)) sigma_x <- default_sigma_x(X, strata)

  n_per_s  <- tabulate(as.integer(strata))
  n1_per_s <- as.integer(tapply(z, strata, sum))
  pi_i     <- (n1_per_s / n_per_s)[as.integer(strata)]

  ## d = (z - pi)' X
  d <- drop(crossprod(z - pi_i, X))

  ## T = d' Sigma_x^{-1} d, computed via the rotation by Sigma_x^{-1/2}
  rot   <- XtX_pseudoinv_sqrt(sigma_x, mat.is.XtX = TRUE)
  d_y   <- drop(d %*% rot)
  Tstat <- sum(d_y^2)

  ## Eigenvalue spectrum for the asymptotic backends:
  ##   lambda = eigenvalues of (Sigma_x^{-1/2} V_d Sigma_x^{-1/2}) = Cov(d_y)
  V_d <- randomization_cov_d(X, strata, n1_per_s)
  V_y <- t(rot) %*% V_d %*% rot
  V_y <- (V_y + t(V_y)) / 2  # symmetrize against numerical drift
  lambda <- Re(eigen(V_y, symmetric = TRUE, only.values = TRUE)$values)
  lambda <- lambda[lambda > sqrt(.Machine$double.eps) * max(lambda, 1)]

  if (null == "simulate") {
    nd <- simulate_T_under_null(X, strata, n1_per_s, sigma_x, n_simulate)
    p  <- (sum(nd >= Tstat) + 1) / (length(nd) + 1)  # plus-one smoothing
    return(list(statistic    = Tstat,
                p.value      = p,
                df_eff       = length(lambda),
                lambda       = lambda,
                sigma_x_used = sigma_x,
                null_draws   = nd))
  }
  if (null == "satterthwaite_finite") {
    mom <- sigma_x_T_moments(X, strata, n1_per_s, sigma_x)
    p   <- sigma_x_pvalue(Tstat, EM = mom$EM, VarM = mom$VarM,
                          method = "satterthwaite_finite")
    df_eff <- 2 * mom$EM^2 / mom$VarM
  } else {
    p <- sigma_x_pvalue(Tstat, lambda = lambda, method = null)
    df_eff <- if (length(lambda) > 0) sum(lambda)^2 / sum(lambda^2) else 0
  }
  list(statistic    = Tstat,
       p.value      = p,
       df_eff       = df_eff,
       lambda       = lambda,
       sigma_x_used = sigma_x)
}

## Draw one treatment vector from within-stratum simple random sampling
## without replacement, with the per-stratum treated counts fixed at
## n1_per_stratum.  Tested directly in test.sigma_x_test.R.
##
## by_s is a precomputed list of row-index vectors per stratum (so the
## caller can hoist the split() call out of a loop).
draw_within_stratum_z <- function(by_s, n1_per_stratum, N) {
  z <- numeric(N)
  for (s in seq_along(by_s)) {
    idx     <- by_s[[s]]
    treated <- sample.int(length(idx), n1_per_stratum[s])
    z[idx[treated]] <- 1
  }
  z
}

## Draw `n_simulate` random within-stratum assignments and compute T for each.
## All draws are uniform over the within-stratum simple random sample without
## replacement (matching SRSWOR), so the empirical distribution converges to
## the exact randomization null at the rate 1/sqrt(B).
simulate_T_under_null <- function(X, strata, n1_per_stratum, sigma_x,
                                   n_simulate) {
  X <- as.matrix(X)
  N <- nrow(X)
  by_s    <- split(seq_len(N), strata)
  n_per_s <- tabulate(as.integer(strata))
  pi_i    <- (n1_per_stratum / n_per_s)[as.integer(strata)]
  rot     <- XtX_pseudoinv_sqrt(sigma_x, mat.is.XtX = TRUE)
  Y       <- as.matrix(X %*% rot)

  out <- numeric(n_simulate)
  for (b in seq_len(n_simulate)) {
    z   <- draw_within_stratum_z(by_s, n1_per_stratum, N)
    d_y <- drop(crossprod(z - pi_i, Y))
    out[b] <- sum(d_y^2)
  }
  out
}

##' Run the sigma_x test on one stratification of an aggregated design
##'
##' Hook used by \code{balanceTest()}.  Extracts the cluster-aggregated raw X,
##' strata factor, and treatment indicator from a \code{DesignOptions}-derived
##' object, drops units with NA strata, and delegates to \code{sigma_x_test()}.
##' @keywords internal
sigma_x_inferentials <- function(aggDesign, stratification_name,
                                  sigma_x = NULL,
                                  null = "satterthwaite_finite",
                                  n_simulate = 1000) {
  X <- as.matrix(aggDesign@Covariates)
  ss_raw <- aggDesign@StrataFrame[, stratification_name]
  z <- as.numeric(aggDesign@Z)
  keep <- !is.na(ss_raw)
  X <- X[keep, , drop = FALSE]
  z <- z[keep]
  strata <- droplevels(factor(ss_raw[keep]))
  sigma_x_test(X, strata, z,
               sigma_x = sigma_x,
               null = null,
               n_simulate = n_simulate)
}
