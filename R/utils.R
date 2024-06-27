## This file contains some small helper functions.

##' Get p-value for Z-stats
##'
##' @param zs A Z-statistic.
##' @return A P-value
makePval <- function(zs) {
  2 * pnorm(abs(zs), lower.tail = FALSE)
}

##' Returns \code{formula} attribute of an \code{xbal} object.
##'
##' @param x An \code{xbal} object.
##' @param ... Ignored.
##' @return The formula corresponding to \code{xbal}.
##' @aliases formula.balancetest
##' @export
formula.xbal <- function(x, ...) {
  attr(x, "fmla")
}

##' Safe way to temporarily override options()
##'
##' @param optionsToChange Which options.
##' @param fun Function to run with new options.
##' @return Result of \code{fun}.
withOptions <- function(optionsToChange, fun) {
  oldOpts <- options()
  options(optionsToChange)
  # store the old values of the options, just for the options that were changed
  oldOptValues <- oldOpts[names(optionsToChange)]
  tryCatch(fun(), finally = options(oldOptValues))
}

## Our own version of these to handle the signif stars.
### print.ftable<-function (x, digits = getOption("digits"), ...) {
###  write.ftable(x, quote = FALSE, digits = digits)
### }
###
### write.ftable<-function (x, file = "", quote = TRUE, append = FALSE, digits = getOption("digits"),justify.labels="right",justify.data="right",...)
### {
###    r <- RItools:::format.ftable(x, quote = quote, digits = digits,justify.labels=justify.labels,justify.data=justify.data,...)
###    cat(t(r), file = file, append = append, sep = c(rep(" ",
###        ncol(r) - 1), "\n"))
###    invisible(x)
### }
###
### format.ftable<-function (x, quote = TRUE, digits = getOption("digits"), justify.labels="left",justify.data="right", ...)
### {
###    if (!inherits(x, "ftable"))
###        stop("'x' must be an \"ftable\" object")
###    charQuote <- function(s) if (quote)
###        paste("\"", s, "\"", sep = "")
###    else s
###    makeLabels <- function(lst) {
###        lens <- sapply(lst, length)
###        cplensU <- c(1, cumprod(lens))
###        cplensD <- rev(c(1, cumprod(rev(lens))))
###        y <- NULL
###        for (i in rev(seq_along(lst))) {
###            ind <- 1 + seq.int(from = 0, to = lens[i] - 1) *
###                cplensD[i + 1]
###            tmp <- character(length = cplensD[i])
###            tmp[ind] <- charQuote(lst[[i]])
###            y <- cbind(rep(tmp, times = cplensU[i]), y)
###        }
###        y
###    }
###    makeNames <- function(x) {
###        nmx <- names(x)
###        if (is.null(nmx))
###            nmx <- rep("", length.out = length(x))
###        nmx
###    }
###    xrv <- attr(x, "row.vars")
###    xcv <- attr(x, "col.vars")
###    LABS <- cbind(rbind(matrix("", nrow = length(xcv), ncol = length(xrv)),
###        charQuote(makeNames(xrv)), makeLabels(xrv)), c(charQuote(makeNames(xcv)),
###        rep("", times = nrow(x) + 1)))
###    DATA <- rbind(if (length(xcv))
###        t(makeLabels(xcv)), rep("", times = ncol(x)), format(unclass(x),
###        digits = digits))
###    cbind(apply(LABS, 2, format, justify = justify.labels), apply(DATA,
###        2, format, justify = justify.data))
### }

##' Select variables, strata, and statistics from a \code{xbal} or \code{balancetest} object
##'
##' If any of the arguments are not specified, all the of relevant
##' items are included.
##'
##' @param x The \code{xbal} object, the result of a call to
##'   \code{\link{xBalance}} or \code{\link{balanceTest}}
##' @param vars The variable names to select.
##' @param strata The strata names to select.
##' @param stats The names of the variable level statistics to select.
##' @param tests The names of the group level tests to select.
##' @param ... Other arguments (ignored)
##'
##' @return A \code{xbal} object with just the appropriate items
##'   selected.
##' @aliases subset.balancetest
##'
##' @export
subset.xbal <- function(x,
                        vars = NULL,
                        strata = NULL,
                        stats = NULL,
                        tests = NULL,
                        ...) {
  res.dmns <- dimnames(x$results)

  if (is.null(strata)) {
    strata <- res.dmns$strata
  }

  if (is.null(vars)) {
    vars <- res.dmns$vars
  }

  if (is.null(stats)) {
    stats <- res.dmns$stat
  }

  if (is.null(tests)) {
    tests <- colnames(x$overall)
  }

  res <- x$results[vars, stats, strata, drop = F]
  ovr <- x$overall[strata, tests, drop = F]

  if (!is.null(ovr)) {
    attr(ovr, "tcov") <- attr(x$overall, "tcov")[strata]
  }

  keep_this_var <- res.dmns$vars %in% vars
  attr(res, "NMpatterns") <- attr(x$res, "NMpatterns")[keep_this_var]
  attr(res, "originals") <- attr(x$results, "originals")[keep_this_var]
  attr(res, "term.labels") <- attr(x$results, "term.labels")
  attr(res, "include.NA.flags") <- attr(x$results, "include.NA.flags")


  tmp <- list(results = res, overall = ovr)
  class(tmp) <- c("xbal", "list")
  attr(tmp, "fmla") <- attr(x, "fmla")
  attr(tmp, "report") <- attr(x, "report")

  return(tmp)
}


## SparseM-related

## Turn a factor variable into a sparse matrix of 0's and 1's, such that if observation i
## has the jth level then there is a 1 at position (i,j) (but nowhere else in row i).
##
## NA's give rise to rows with no 1s.
## As the result is only meaningful in the context of the SparseM package,
## function requires that SparseM be loaded.
## @title Sparse matrix dummy coding of a factor variable (omitting the intercept)
## @param thefactor Factor variable, or object inheriting from class factor
## @return Sparse csr matrix the columns of which are dummy variables for levels of thefactor
## @author Ben Hansen
## @examples
## sparse_mod_matrix <-  SparseMMFromFactor(iris$Species)
## mod_matrix <- model.matrix(~Species-1, iris)
## all.equal(as.matrix(sparse_mod_matrix),
##           mod_matrix, check.attributes=FALSE)
SparseMMFromFactor <- function(thefactor) {
  stopifnot(inherits(thefactor, "factor"))
  theNA <- ## if (inherits(thefactor, "optmatch")) !matched(thefactor) else
    is.na(thefactor)

  if (all(theNA)) stop("No non-NA's in thefactor")

  nlev <- nlevels(thefactor)
  nobs <- length(thefactor)
  theint <- as.integer(thefactor)
  if (any(theNA)) theint[theNA] <- 1L # odd; but note how we compensate in def of ra slot below
  new("matrix.csr",
    ja = theint,
    ia = 1L:(nobs + 1L),
    ra = (1L - theNA),
    dimension = c(nobs, nlev) #+sum(theNA)
  )
}

#' SparseM::slm.fit.csr, made tolerant to faults that recur in RItools
#' 
#' [SparseM's slm.fit.csr()] expects a full-rank x that's not just
#' a column of 1s. This variant somewhat relaxes these expectations.
#'
#' `slm.fit.csr` has a bug for intercept only models
#' (admittedly, these are generally a little silly to be done as a
#' sparse matrix), but in order to avoid duplicate code, if
#' everything is in a single strata, we use the intercept only model.
#' 
#' This function's expectation of x is that either it has full column
#' rank, or the reduced submatrix of x that excludes all-zero columns
#' has full column rank. (When this expectation is not met, it's
#' likely that [SparseM::chol()] will fail, causing this function to
#' error; the error messages won't necessarily suggest this.) The
#' positions of nonzero x-columns (ie columns with nonzero entries)
#' are returns as the value of `gramian_reduction_index`, while `chol`
#' is the Cholesky decomposition of that submatrix's Gramian.
#' 
#' @param x As slm.fit.csr
#' @param y As slm.fit.csr
#' @param ... As slm.fit.csr
#' @return A list consisting of:
#'   \item{coefficients}{coefficients}
#'   \item{chol}{Cholesky factor of Gramian matrix \eqn{x'x}}
#'   \item{residuals}{residuals}
#'   \item{fitted}{fitted values}
#'   \item{df.residual}{degrees of freedom}
#'   \item{gramian_reduction_index}{Column indices identifying reduction of x matrix of which Gramian is taken; see Details}

slm_fit_csr <- function(x, y, ...) {
  if (is.matrix(y)) {
    n <- nrow(y)
    ycol <- ncol(y)
  } else {
    n <- length(y)
    ycol <- 1
  }
  p <- x@dimension[2]
  if (n != x@dimension[1]) {
    stop("x and y don't match n")
  }

  temp_sol <- SparseM_solve(x, y, ...)
  coef <- temp_sol[["coef"]]
  chol <- temp_sol[["chol"]]

  if (is.vector(coef)) {
    coef <- matrix(coef, ncol = ycol, nrow = p)
  }
  fitted <- as.matrix(x %*% coef)
  resid <- y - fitted
  df <- n - p
  list(
    coefficients = coef,
    chol = chol,
    residuals = resid,
    fitted = fitted, 
    df.residual = df,
    gramian_reduction_index = temp_sol[["gramian_reduction_index"]]
  )
}


#' Helper function to slm_fit_csr
#' 
#' This function generates a matrix that can be used to reduce
#' the dimensions of x'x and xy such that positive definiteness is
#' ensured and more practically, that SparseM::chol will work
#' 
#' @param zeroes logical vector indicating which entries of the diagonal of x'x are zeroes.
#' @return SparseM matrix that will reduce the dimension of x'x and xy 
#' @importFrom SparseM chol backsolve
gramian_reduction <- function(zeroes)
{
  if (all(zeroes))
  {
    stop("Diagonal of X'X is all zeroes. Unable to proceed.")
  }
  
  num_rows <- length(zeroes)
  num_cols <- sum(!zeroes)
  non_zero_indices <- which(!zeroes)
  
  # Calculate the column indices for non-zero values
  col_indices <- sapply(non_zero_indices, 
                        function(i) i - sum(zeroes[1:i]))
  values <- rep(1, num_cols)
  
  # Define the row pointer array 
  ia <- cumsum(c(1, !zeroes))
  
  dimension <- as.integer(c(num_rows, num_cols))
  reducing_matrix <- new("matrix.csr", ra = values, 
                         ja = as.integer(col_indices), 
                         ia = as.integer(ia), 
                         dimension = dimension)
  
  return(reducing_matrix)
}


#' Helper function to slm_fit_csr
#' 
#' This function performs some checks and takes action to 
#' ensure positive definiteness of matrices passed to SparseM functions.
#' 
#' @param x A slm.fit.csr
#' @param y A slm.fit.csr
#' @param ... A slm.fit.csr
#' @return list containing coefficients (vector or matrix), the Cholesky decomposition (of class matrix.csr.chol), and a vector specifying the indices of which values on the diagonal of x'x are nonzero. These are named "coef", "chol" and "gramian_reduction_index", respectively.
#' @importFrom SparseM chol backsolve
SparseM_solve <- function(x, y, ...)
{
  xy <- t(x) %*% y
  xprimex <- t(x) %*% x
  diag.xx <- diag(xprimex)
  zeroes <- diag.xx == 0
  if (any(zeroes)) #check explicitly for zeroes here so we don't do matrix math without needing to
  { # this branch deals with issue 134
    reducing_matrix <- gramian_reduction(zeroes)
    xpx.sub <- t(reducing_matrix) %*% xprimex %*% reducing_matrix
    xy.sub <- t(reducing_matrix) %*% xy
    chol.result <- SparseM::chol(xpx.sub, ...)
    coef.nonzero <- SparseM::backsolve(chol.result, xy.sub)
    num_rows <- length(zeroes)
    coef.all <- numeric(num_rows)
    coef.all[!zeroes] <- coef.nonzero
  } else
  {
    chol.result <- SparseM::chol(xprimex, ...)
    coef.all <- SparseM::backsolve(chol.result, xy)
  }
  
  return(list("coef" = coef.all, 
              "chol" = chol.result, 
              "gramian_reduction_index" = which(!zeroes)))
}

## slm.wfit with two fixes
##
## slm.wfit shares the intercept-only issue with slm.fit,
## and in addition has an issue where it carries forward
## incorrect residuals and fitted values.
##
## @param x As slm.wfit
## @param y As slm.wfit
## @param weights As slm.wfit
## @param ... As slm.wfit
## @return As slm.wfit
##' @importFrom SparseM is.matrix.csr
slm.wfit.csr <- function(x, y, weights, ...) {
  if (!is.matrix.csr(x)) {
    stop("model matrix must be in sparse csr mode")
  }
  if (!is.numeric(y)) {
    stop("response must be numeric")
  }
  if (any(weights < 0)) {
    stop("negative weights not allowed")
  }
  contr <- attr(x, "contrasts")
  w <- sqrt(weights)
  wx <- as(w, "matrix.diag.csr") %*% x
  wy <- y * w
  fit <- slm_fit_csr(wx, wy, ...)

  fit$fitted <- as.matrix(x %*% fit$coef)
  fit$residuals <- y - fit$fitted

  fit$contrasts <- attr(x, "contrasts")
  fit
}

### Other linear algebra
## Modeled on MASS's \code{ginv}
##
##
## @title Matrix square root of XtX's pseudoinverse
## @param mat double-precision matrix
## @param mat.is.XtX is mat a crossproduct of an X matrix, or X itself?
## @param tol tolerance
## @return matrix of \code{ncol(mat)} rows and col rank (mat) columns
## @author Ben Hansen
## @keywords internal
XtX_pseudoinv_sqrt <- function(mat, mat.is.XtX = FALSE, tol = .Machine$double.eps^0.5) {
  pst.svd <- try(svd(mat, nu = 0))

  if (inherits(pst.svd, "try-error")) {
    pst.svd <- propack.svd(mat)
  }

  d <- if (mat.is.XtX) sqrt(pst.svd$d) else pst.svd$d
  v <- pst.svd$v

  Positive <- d > max(tol * d[1], 0)
  Positive[is.na(Positive)] <- FALSE

  if (all(Positive)) {
    ytl <- v *
      matrix(1 / d, nrow = ncol(mat), ncol = length(d), byrow = T)
  } else if (!any(Positive)) {
    ytl <- array(0, c(ncol(mat), 0))
  } else {
    ytl <- v[, Positive, drop = FALSE] *
      matrix(1 / d[Positive], ncol = sum(Positive), nrow = ncol(mat), byrow = TRUE)
  }
  r <- sum(Positive)
  attr(ytl, "r") <- r

  return(ytl)
}
