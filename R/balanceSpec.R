##' Specify a balance test against a supplied covariance matrix
##'
##' A \code{balanceSpec} packages what \code{\link{balanceTest}} needs to test
##' covariate balance when the covariance of the covariates is supplied rather
##' than calculated from the data at hand: the treatment assignment, the matched
##' sets or other strata, an optional clustering variable and unit weights, the
##' covariance matrix, and the variables to test. Passing the object to
##' \code{balanceTest()} along with a data frame runs the test.
##'
##' The supplied matrix stands in for the within-stratum covariance of the
##' covariates, taken to be the same in every stratum. Under that reading the
##' null covariance of the vector of adjusted treatment-control differences is
##' the supplied matrix times a scalar that depends only on the design (the
##' stratum sizes, their treated and control counts, and the stratum weights),
##' so the inverse square root the omnibus test needs is computed once, when the
##' specification is made, and reused by every test run against it.
##'
##' The columns of \code{covariance} must be named for the expanded covariate
##' columns, as \code{\link{expandCovariates}} produces them: factors contribute
##' one column per level. A test stops with a message naming any expanded
##' covariate column the covariance lacks. Covariance columns that the data do
##' not produce (for instance a factor level absent from this sample) are set
##' aside, and the inverse is recomputed for the remaining block.
##'
##' Missing covariate values are imputed to the stratum mean, as in the formula
##' method, and the null covariance is reduced accordingly: each entry of the
##' supplied matrix is scaled by the share of units, stratum by stratum, with
##' both of its variables observed. With any missing values that scaling is not
##' a scalar, so the inverse is recomputed for the test at hand rather than
##' reused. Unit weights other than 1 and clusters are passed through to the
##' design, but the design scalar does not account for them: in that case the
##' supplied covariance has to describe the weighted (or cluster-aggregated)
##' rows.
##'
##' The not-missing indicator rows that \code{include.NA.flags} adds are
##' reported descriptively but not tested, and do not enter the omnibus
##' statistic: the supplied covariance has no entries for them. A unit with
##' zero weight counts as unobserved on every variable, which is what its
##' zero row amounts to in \code{HB08}'s sum of products.
##'
##' \code{subset} and \code{unit.weights} are not arguments of this method:
##' subset the data frame beforehand, and give unit weights in the
##' specification.
##'
##' @param treatment The treatment assignment: the name of a column of the data
##'   frame the test will be run on, or a vector with one entry per row.
##' @param covariance A symmetric numeric matrix with matching row and column
##'   names, the covariance of the expanded covariates.
##' @param strata The matched sets or strata: a column name, a vector of column
##'   names (one stratification each), a factor or vector with one entry per
##'   row, or \code{NULL} for no stratification.
##' @param variables The variables to test: a character vector of column names
##'   or a one-sided formula, in terms of the data frame's columns (not the
##'   expanded columns). \code{NULL}, the default, takes every column not used
##'   for treatment, strata, cluster or weights.
##' @param cluster Optional clustering variable, a column name or a vector.
##' @param unit.weights Optional unit weights, a column name or a numeric vector.
##' @param inverse_sqrt Optional precomputed inverse square root of
##'   \code{covariance}, as \code{XtX_pseudoinv_sqrt(covariance, mat.is.XtX =
##'   TRUE)} returns it, with its \code{r} attribute (the rank). Computed here
##'   when not supplied.
##' @param tol Relative tolerance for the pseudo-inverse: singular values of
##'   the covariance below \code{tol} times the largest count as zero. (This
##'   is applied to the covariance itself, not to the square roots of its
##'   singular values as \code{XtX_pseudoinv_sqrt()}'s own \code{tol} is,
##'   so the rank of a singular matrix, such as one over a factor's full set
##'   of indicator columns, does not depend on rounding at machine precision.)
##' @return An object of class \code{balanceSpec}.
##' @seealso \code{\link{balanceTest}}, \code{\link{expandCovariates}}
##' @export
##' @examples
##' data(nuclearplants)
##' Sigma <- cov(expandCovariates(~ date + t1 + t2 + cap, nuclearplants))
##' spec <- balanceSpec(treatment = "pr", strata = "pt", covariance = Sigma,
##'                     variables = c("date", "t1", "t2", "cap"))
##' balanceTest(spec, data = nuclearplants)
balanceSpec <- function(treatment, covariance, strata = NULL, variables = NULL,
                        cluster = NULL, unit.weights = NULL,
                        inverse_sqrt = NULL, tol = .Machine$double.eps^0.5) {
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance)) {
    stop("covariance must be a square numeric matrix")
  }
  if (is.null(colnames(covariance)) ||
      !identical(rownames(covariance), colnames(covariance))) {
    stop("covariance must have matching row and column names")
  }
  if (anyDuplicated(colnames(covariance))) {
    stop("covariance has repeated column names")
  }
  if (anyNA(covariance)) {
    stop("covariance has missing entries")
  }
  if (!isTRUE(all.equal(covariance, t(covariance), check.attributes = FALSE))) {
    stop("covariance must be symmetric")
  }
  if (!is.null(variables) && !is.character(variables) && !inherits(variables, "formula")) {
    stop("variables must be a character vector, a one-sided formula, or NULL")
  }
  if (is.null(inverse_sqrt)) {
    inverse_sqrt <- .covarianceInverseSqrt(covariance, tol)
  } else {
    r <- attr(inverse_sqrt, "r")
    if (!is.matrix(inverse_sqrt) || nrow(inverse_sqrt) != ncol(covariance) ||
        is.null(r) || ncol(inverse_sqrt) != r) {
      stop("inverse_sqrt must have one row per column of covariance and as many ",
           "columns as its 'r' attribute, as XtX_pseudoinv_sqrt() returns it")
    }
    if (!is.null(rownames(inverse_sqrt)) &&
        !identical(rownames(inverse_sqrt), colnames(covariance))) {
      stop("inverse_sqrt rows are not in the order of the covariance's columns")
    }
    .covarianceInverseSqrt(covariance, tol, check.only = TRUE)
  }
  structure(list(treatment = treatment,
                 strata = strata,
                 cluster = cluster,
                 unit.weights = unit.weights,
                 variables = variables,
                 covariance = covariance,
                 inverse_sqrt = inverse_sqrt,
                 tol = tol),
            class = "balanceSpec")
}

##' Pseudo-inverse square root of a covariance, from its eigendecomposition
##'
##' One symmetric eigendecomposition serves both the positive semi-definite
##' check and the inverse. Eigenvalues below \code{tol} times the largest
##' count as zero, a cut on the matrix's own scale; the square roots that
##' \code{XtX_pseudoinv_sqrt(mat, mat.is.XtX = TRUE)} compares against its
##' tolerance would put that cut at machine precision and let the rank of a
##' singular matrix vary with the BLAS.
##'
##' @param covariance A symmetric matrix.
##' @param tol Relative tolerance on the eigenvalues of \code{covariance}.
##' @param check.only Only check positive semi-definiteness, invisibly
##'   returning \code{NULL}?
##' @return A matrix \code{V} with as many rows as \code{covariance} and one
##'   column per positive eigenvalue, such that \code{V \%*\% t(V)} is the
##'   pseudo-inverse, with the rank as its \code{r} attribute, in the form
##'   \code{XtX_pseudoinv_sqrt()} returns.
##' @keywords internal
.covarianceInverseSqrt <- function(covariance, tol, check.only = FALSE) {
  e <- eigen(covariance, symmetric = TRUE, only.values = check.only)
  lambda <- e$values
  top <- if (length(lambda)) max(abs(lambda), 0) else 0
  if (length(lambda) && min(lambda) < -tol * max(top, 1)) {
    stop("covariance must be positive semi-definite (smallest eigenvalue ",
         format(min(lambda)), ")")
  }
  if (check.only) return(invisible(NULL))
  positive <- lambda > tol * top
  positive[is.na(positive)] <- FALSE
  ans <- e$vectors[, positive, drop = FALSE] *
    rep(1 / sqrt(lambda[positive]), each = nrow(covariance))
  attr(ans, "r") <- sum(positive)
  ans
}

##' @method print balanceSpec
##' @export
print.balanceSpec <- function(x, ...) {
  describe <- function(v) {
    if (is.null(v)) return("none")
    if (is.character(v) && length(v) <= 3L && all(nchar(v) < 40)) {
      return(paste(sprintf("'%s'", v), collapse = ", "))
    }
    sprintf("<vector of length %d>", length(v))
  }
  vars <- if (is.null(x$variables)) {
    "all columns not otherwise used"
  } else if (inherits(x$variables, "formula")) {
    paste(deparse(x$variables), collapse = " ")
  } else {
    paste(x$variables, collapse = ", ")
  }
  cat("Balance test specification with a supplied covariance\n")
  cat("  treatment:   ", describe(x$treatment), "\n", sep = "")
  cat("  strata:      ", describe(x$strata), "\n", sep = "")
  cat("  cluster:     ", describe(x$cluster), "\n", sep = "")
  cat("  unit.weights:", describe(x$unit.weights), "\n", sep = "")
  cat("  variables:   ", vars, "\n", sep = "")
  cat(sprintf("  covariance:  %d by %d, rank %d\n",
              nrow(x$covariance), ncol(x$covariance), attr(x$inverse_sqrt, "r")))
  invisible(x)
}

##' @rdname balanceTest
##' @order 3
##' @param unstratified For the \code{balanceSpec} method: also report the
##'   unstratified comparison (the \code{--} stratification)?
##' @method balanceTest balanceSpec
##' @export
balanceTest.balanceSpec <- function(fmla, data,
                                    unstratified = TRUE,
                                    include.NA.flags = TRUE,
                                    p.adjust.method = "holm",
                                    stratum.weights = harmonic_times_mean_weight,
                                    covariate.scales = NULL,
                                    post.alignment.transform = NULL,
                                    ...) {
  .stopOnUnusedArguments(match.call(expand.dots = FALSE)$...)
  spec <- fmla
  if (!is.data.frame(data)) stop("data must be a data frame")
  stopifnot(is.null(post.alignment.transform) || is.function(post.alignment.transform))
  n <- nrow(data)

  ## Each structural piece is a column name (or several, for strata) or a
  ## vector to be added as a column.
  resolve <- function(value, what, newname) {
    if (is.null(value)) return(character(0))
    if (is.character(value) && length(value) < n && all(value %in% colnames(data))) return(value)
    if (length(value) != n) {
      stop(sprintf(paste("%s is neither column name(s) of data (%s)",
                         "nor a vector with one entry per row (%d rows)"),
                   what, paste(value[seq_len(min(3L, length(value)))], collapse = ", "), n))
    }
    data[[newname]] <<- value
    newname
  }
  zname <- resolve(spec$treatment, "treatment", ".treatment")
  if (length(zname) != 1L) stop("treatment must name one column or be one vector")
  snames <- resolve(spec$strata, "strata", ".strata")
  cname <- resolve(spec$cluster, "cluster", ".cluster")
  if (length(cname) > 1L) stop("at most one cluster variable is allowed")
  wname <- resolve(spec$unit.weights, "unit.weights", ".unit.weights")
  if (length(wname) > 1L) stop("unit.weights must name one column or be one vector")

  if (!unstratified && !length(snames)) {
    stop("unstratified = FALSE needs strata in the specification")
  }
  structural <- c(zname, snames, cname, wname)
  rhs <- if (is.null(spec$variables)) {
    vars <- setdiff(colnames(data), structural)
    if (!length(vars)) stop("no variables left to test")
    paste(sprintf("`%s`", vars), collapse = " + ")
  } else if (inherits(spec$variables, "formula")) {
    paste(deparse(spec$variables[[length(spec$variables)]]), collapse = " ")
  } else {
    missing_vars <- setdiff(spec$variables, colnames(data))
    if (length(missing_vars)) {
      stop("variables not in data: ", paste(missing_vars, collapse = ", "))
    }
    paste(sprintf("`%s`", spec$variables), collapse = " + ")
  }
  fmla <- as.formula(paste0(
    "`", zname, "` ~ ", rhs,
    if (length(snames)) paste0(" + strata(`", snames, "`)", collapse = "") else "",
    if (length(cname)) paste0(" + cluster(`", cname, "`)") else "",
    if (unstratified) "" else " - 1"),
    env = parent.frame())

  w <- if (length(wname)) data[[wname]] else rep(1, n)
  if (!is.numeric(w)) w <- as.numeric(w)
  data$'(weights)' <- w

  if (is.null(covariate.scales)) {
    ## standard differences over the supplied dispersion, where it is usable;
    ## a zero or non-finite entry leaves that variable to the sample estimate
    scales <- diag(spec$covariance)
    names(scales) <- colnames(spec$covariance)
    scales <- scales[is.finite(scales) & scales > 0]
    covariate.scales <- sqrt(scales)
  }

  ans <- .balanceTestFromFrame(fmla, data,
                               stratum.weights = stratum.weights,
                               include.NA.flags = include.NA.flags,
                               covariate.scales = covariate.scales,
                               post.alignment.transform = post.alignment.transform,
                               inferentials.calculator = HB08_fixed(spec),
                               p.adjust.method = p.adjust.method)
  attr(ans, "covariance") <- "supplied"
  ans
}

##' Hansen and Bowers (2008) inferentials against a supplied covariance
##'
##' Returns a function with the interface of \code{\link{HB08}}, for use as
##' \code{balanceTest}'s \code{inferentials.calculator}, that takes the null
##' covariance of the adjusted differences to be the specification's
##' covariance matrix times a scalar fixed by the design: the sum over strata
##' of the squared stratum weight ratio times the reciprocal of
##' \eqn{1/n_{tb} + 1/n_{cb}}, the stratum's treated and control counts.
##' Strata with one unit, or lacking either group, contribute nothing, as in
##' \code{HB08}. Where covariate values were imputed, the stratum's
##' contribution to an entry is further scaled by the share of its units with
##' both variables observed, which is what the sample covariance of the
##' imputed columns in \code{HB08} amounts to.
##'
##' The returned function keeps state: the message about covariance columns
##' the data do not produce is given once, and the inverse for a reduced
##' block is kept for later stratifications with the same block.
##'
##' @param spec A \code{\link{balanceSpec}}.
##' @return A function of a \code{CovsAlignedToADesign} returning a list as
##'   \code{\link{HB08}} does.
##' @keywords internal
HB08_fixed <- function(spec) {
  Sigma <- spec$covariance
  inv_full <- spec$inverse_sqrt
  tol <- spec$tol

  cache <- list(tested = NULL, inv = NULL, told = FALSE)

  function(alignedcovs) {
    zz <- as.numeric(alignedcovs@Z)
    S <- alignedcovs@StrataMatrix
    Covs <- alignedcovs@Covariates
    NMcov <- alignedcovs@NotMissing
    n_ <- nrow(Covs)
    wtr <- alignedcovs@StrataWeightRatio
    stopifnot(length(wtr) == n_, identical(dim(NMcov), dim(Covs)))
    vars <- colnames(Covs)

    ## which columns the supplied covariance covers
    covariate_cols <- vars[!startsWith(vars, "(")]
    not_in_sigma <- setdiff(covariate_cols, colnames(Sigma))
    if (length(not_in_sigma)) {
      stop("expanded covariate column(s) not in the supplied covariance: ",
           paste(not_in_sigma, collapse = ", "),
           ". The covariance must be over the columns expandCovariates() produces.")
    }
    tested <- intersect(colnames(Sigma), vars)
    if (!length(tested)) {
      stop("none of the supplied covariance's columns is among the covariates")
    }
    if (length(tested) < ncol(Sigma) && !cache$told) {
      message(sprintf(paste("%d column(s) of the supplied covariance are not among the",
                            "covariates (%s); testing the remaining %d and recomputing",
                            "the inverse for them"),
                      ncol(Sigma) - length(tested),
                      paste(setdiff(colnames(Sigma), tested), collapse = ", "),
                      length(tested)))
      cache$told <<- TRUE
    }

    ## stratum counts and the design scalar, stratum by stratum
    n <- sparseToVec(t(S) %*% rep(1, n_))
    n1 <- sparseToVec(t(S) %*% zz)
    n0 <- n - n1
    n1_over_n <- sparseToVec(S %*% (n1 / n))
    wtr_b <- sparseToVec(t(S) %*% wtr) / n
    usable <- n > 1 & n0 > 0 & n1 > 0
    c_b <- ifelse(usable, wtr_b^2 / (1 / n0 + 1 / n1), 0)
    cc <- sum(c_b)

    ## adjusted differences, as in HB08
    x_tilde <- Covs * wtr
    ssn <- drop(crossprod(zz - n1_over_n, x_tilde))
    names(ssn) <- vars

    ## Null covariance: sum over strata of c_b times the supplied covariance,
    ## each entry scaled by the stratum's share of units observed on both
    ## variables. With nothing missing every share is 1 and this is cc * Sigma.
    complete <- all(NMcov[, tested] == 1)
    if (complete) {
      scale_mat <- matrix(cc, length(tested), length(tested))
    } else {
      c_row <- sparseToVec(S %*% (c_b / n))
      scale_mat <- crossprod(NMcov[, tested, drop = FALSE] * sqrt(c_row))
    }
    tcov <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
    tcov[tested, tested] <- scale_mat * Sigma[tested, tested]
    ssvar <- diag(tcov)
    zero_variance <- ssvar <= .Machine$double.eps
    zstat <- ifelse(zero_variance, NA_real_, ssn / sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

    ## the inverse square root: the precomputed one whenever it applies
    if (complete && length(tested) == ncol(Sigma)) {
      inv <- inv_full
      divisor <- sqrt(cc)
    } else if (complete && identical(tested, cache$tested)) {
      inv <- cache$inv
      divisor <- sqrt(cc)
    } else if (complete) {
      inv <- .covarianceInverseSqrt(Sigma[tested, tested, drop = FALSE], tol)
      cache$tested <<- tested
      cache$inv <<- inv
      divisor <- sqrt(cc)
    } else {
      inv <- .covarianceInverseSqrt(tcov[tested, tested, drop = FALSE], tol)
      divisor <- 1
    }

    check_for_degenerate(inv, n_, ncol(S))

    DF <- attr(inv, "r")
    csq <- if (cc > 0 && DF > 0) {
      mvz <- drop(crossprod(ssn[tested], inv)) / divisor
      sum(mvz^2)
    } else NA_real_

    list(z = zstat, p = p, Msq = csq, DF = DF,
         adj.diff.of.totals = ssn, tcov = tcov)
  }
}
