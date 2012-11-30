##This file contains some small helper functions.

makePval<-function(zs){2*pnorm(abs(zs),lower=FALSE)}

formula.xbal<-function(x,...){
  attr(x,"fmla")
}

# ============================================================================
# = withOptions helper provides a safe way to temporarily override options() =
# ============================================================================
withOptions <- function(optionsToChange, fun) {
  oldOpts <- options()
  do.call(options, optionsToChange)
  tryCatch(fun(), finally = do.call(options, oldOpts))
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

#' Select variables, strata, groups, and statistics from a \code{xbal} object
#'
#' If any of the arguments are not specified, all the of relevant items are
#' included.
#'
#' @param xbal The \code{xbal} object, the result of a call to
#' \code{\link{xBalance}}
#' @param which.vars The variable names to select.
#' @param which.strata The strata names to select.
#' @param which.stats The names of the variable level statistics to select.
#' @param which.groups The names of the groups to select.
#' @param which.tests The names of the group level tests to select.
#'
#' @return A \code{xbal} object with just the appropriate items selected.
#' @export
select <- function(xbal, 
                   which.vars   = NULL, 
                   which.strata = NULL, 
                   which.stats  = NULL, 
                   which.groups = NULL,
                   which.tests  = NULL) {

  res.dmns <- dimnames(xbal$results)
  grp.dmns <- dimnames(xbal$groups)

  if (is.null(which.strata)) {
    which.strata <- res.dmns$strata
  }

  if (is.null(which.vars)) {
    which.vars <- res.dmns$vars
  }

  if (is.null(which.stats)) {
    which.stats <- res.dmns$stat
  }

  if (is.null(which.groups)) {
    which.groups <- grp.dmns$groups
  }

  if (is.null(which.tests)) {
    which.tests <- grp.dmns$tests
  }

  res <- xbal$results[which.vars, which.stats, which.strata, drop = F] 
  grp <- xbal$groups[which.strata, which.tests, , drop = F] 

  return(list(results = res, groups = grp))  
}
