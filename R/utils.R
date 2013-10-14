##This file contains some small helper functions.

makePval<-function(zs){2*pnorm(abs(zs),lower.tail=FALSE)}

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

#' Fill in an array from a function that takes the indices as arguments.
#'
#' @param f The function to be called.
#' @param p The parameter list(a = c(1,2,3), b = c(7,8), c = c(9,10,11)) will create a 3x2x3 array.
#' @return An array of the cartesian product of the parameters, with f applied to each combination.
farray <- function(f, p) {

  p <- lapply(p, sort)
  n <- length(p)
  j <- sapply(p, length)
  k <- prod(j)
  
  args <- make_args_mtx(p)
  x <- vapply(args, FUN.VALUE = numeric(1), function(a) do.call(f, a))

  array(x, dim = j, dimnames = p)
}

### From:
### http://stackoverflow.com/questions/6192848/how-to-generalize-outer-to-n-dimensions
list_args <- Vectorize(function(a,b) c(as.list(a), as.list(b)), SIMPLIFY = FALSE)

make_args_mtx <- function(alist) {
  Reduce(function(x, y) outer(x, y, list_args), alist)
}

