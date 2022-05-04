##This file contains some small helper functions.

##' Get p-value for Z-stats
##'
##' @param zs A Z-statistic.
##' @return A P-value
makePval<-function(zs) {
  2*pnorm(abs(zs),lower.tail=FALSE)
}

##' Returns \code{formula} attribute of an \code{xbal} object.
##'
##' @param x An \code{xbal} object.
##' @param ... Ignored.
##' @return The formula corresponding to \code{xbal}.
##' @aliases formula.balancetest
##' @export
formula.xbal<-function(x,...){
  attr(x,"fmla")
}

##' Safe way to temporarily override options()
##'
##' @param optionsToChange Which options.
##' @param fun Function to run with new options.
##' @return Result of \code{fun}.
withOptions <- function(optionsToChange, fun) {
  oldOpts <- options()
  on.exit(options(oldOpts))
  options(optionsToChange)
  tryCatch(fun(), finally = options(oldOpts))
}

##Our own version of these to handle the signif stars.
###print.ftable<-function (x, digits = getOption("digits"), ...) {
###  write.ftable(x, quote = FALSE, digits = digits)
###}
###
###write.ftable<-function (x, file = "", quote = TRUE, append = FALSE, digits = getOption("digits"),justify.labels="right",justify.data="right",...)
###{
###    r <- RItools:::format.ftable(x, quote = quote, digits = digits,justify.labels=justify.labels,justify.data=justify.data,...)
###    cat(t(r), file = file, append = append, sep = c(rep(" ",
###        ncol(r) - 1), "\n"))
###    invisible(x)
###}
###
###format.ftable<-function (x, quote = TRUE, digits = getOption("digits"), justify.labels="left",justify.data="right", ...)
###{
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
###}

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
                        vars   = NULL,
                        strata = NULL,
                        stats  = NULL,
                        tests  = NULL,
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

  if (!is.null(ovr))
      attr(ovr, "tcov") <- attr(x$overall, "tcov")[strata]

  keep_this_var <- res.dmns$vars %in% vars
  attr(res, "NMpatterns") <- attr(x$res, "NMpatterns")[keep_this_var]
  attr(res, "originals") <- attr(x$results, "originals")[keep_this_var]
  attr(res, "term.labels") <- attr(x$results, "term.labels")
  attr(res, "include.NA.flags") <-  attr(x$results, "include.NA.flags")

  
  tmp <- list(results = res, overall = ovr)
  class(tmp) <- c("xbal", "list")
  attr(tmp, "fmla") <- attr(x, "fmla")
  attr(tmp, "report") <-  attr(x, "report") 

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
  theNA <- ##if (inherits(thefactor, "optmatch")) !matched(thefactor) else
    is.na(thefactor)

  if (all(theNA)) stop("No non-NA's in thefactor") 

  nlev <- nlevels(thefactor)
  nobs <- length(thefactor)
  theint <- as.integer(thefactor)
  if (any(theNA)) theint[theNA] <- 1L# odd; but note how we compensate in def of ra slot below
  new("matrix.csr",
      ja=theint,
      ia=1L:(nobs+1L),
      ra=(1L-theNA),
      dimension = c(nobs, nlev) #+sum(theNA)
      )
}


## slm.fit.csr with a fix
##
## SparseM's slm.fit.csr has a bug for intercept only models
## (admittedly, these are generally a little silly to be done as a
## sparse matrix), but in order to avoid duplicate code, if
## everything is in a single strata, we use the intercept only model.
## @param x As slm.fit.csr
## @param y As slm.fit.csr
## @param ... As slm.fit.csr
## @return As slm.fit.csr
##' @importFrom SparseM chol backsolve
slm.fit.csr.fixed <- function (x, y, ...)
{
    if (is.matrix(y))
        {
            n <- nrow(y)
            ycol <- ncol(y)
        } else {
            n <- length(y)
            ycol <- 1
        }
    p <- x@dimension[2]
    if (n != x@dimension[1])
        stop("x and y don't match n")
    chol <- SparseM::chol(t(x) %*% x, ...)
    xy <- t(x) %*% y
    coef <- SparseM::backsolve(chol, xy)

    if (is.vector(coef)) {
      coef <- matrix(coef, ncol = ycol, nrow = p)
  }

    fitted <- as.matrix(x %*% coef)
    resid <- y - fitted
    df <- n - p
    list(coefficients = coef, chol = chol, residuals = resid,
        fitted = fitted, df.residual = df)
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
slm.wfit.csr <- function (x, y, weights, ...) 
{

    if (!is.matrix.csr(x)) 
        stop("model matrix must be in sparse csr mode")
    if (!is.numeric(y)) 
        stop("response must be numeric")
    if (any(weights < 0)) 
        stop("negative weights not allowed")
    contr <- attr(x, "contrasts")
    w <- sqrt(weights)
    wx <- as(w, "matrix.diag.csr") %*% x
    wy <- y * w
    fit <- slm.fit.csr.fixed(wx, wy, ...)

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
XtX_pseudoinv_sqrt <- function(mat, mat.is.XtX = FALSE, tol = .Machine$double.eps^0.5)
{
    pst.svd <- try(svd(mat, nu=0))

  if (inherits(pst.svd, 'try-error')) {
    pst.svd <- propack.svd(mat)
  }

    d  <-  if (mat.is.XtX) sqrt(pst.svd$d) else pst.svd$d
    v  <- pst.svd$v

    Positive <- d > max(tol * d[1], 0)
    Positive[is.na(Positive)] <- FALSE

  if (all(Positive)) { 
    ytl <- v *
      matrix(1/d, nrow = ncol(mat), ncol = length(d), byrow = T)
  } else if (!any(Positive)) {
    ytl <- array(0, c(ncol(mat), 0) )
  } else  {
    ytl <- v[, Positive, drop = FALSE] *
      matrix(1/d[Positive], ncol = sum(Positive), nrow = ncol(mat), byrow = TRUE)
  }
  r <- sum(Positive)
  attr(ytl, "r") = r

  return(ytl)
}
