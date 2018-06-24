#' Tidy methods for balanceTestI()/xBalance() results
#'
#' [RItools::tidy.xbal()] has a parameter `varnames_crosswalk` not shared with 
#' [RItools::glance.xbal()]. It should be a named character vector, the elements
#' of which give names of columns to be returned and the names of which correspond
#' to columns of xbal objects' \sQuote{results} entry.  Its ordering dictates the order
#' of the result.
#' Its default value translates between conventional xbal
#' column names and broom package conventional names:
#'     \item{vars}{variable name}
#'     \item{Control}{mean of LHS variable = 0 group}
#'     \item{Treatment}{ mean of LHS variable = 1 group}
#'     \item{adj.diff}{T - C diff w/ direct standardization for strata if applicable}
#'     \item{std.diff}{adj.diff/pooled.sd}
#'     \item{pooled.sd}{pooled SD}
#'     \item{statistic}{`z` column from the xbal object}
#'     \item{p.value}{`p` column from the xbal object
#' Additional parameters beyond those listed here are ignored (at this time). 
#' @param x object of class `"xbal"`, result of [RItools::balancetest()]
#'    or [RItools::xBalance()]
#' @param strata which stratification to return info about? defaults
#'    to last one specified in originating function call
#' @param varnames_crosswalk character vector of new names for xbal columns, named by the xbal column
#' @return data frame composed of: for `[RItools::tidy()]`, a column of variable labels (`vars`) and 
#'         additional columns of balance-related stats; for `[RItools::glance()]`, scalars describing 
#'         a combined differences test, if found, and otherwise `NULL`.
#' @S3method
tidy.xbal <- function(x, strata=dimnames(x[['results']])[['strata']][1],
                      varnames_crosswalk=c("z"="statistic", "p"="p.value"),
                      ...
                      )
{
    ans <- x[['results']][,,strata, drop=FALSE]
    dim(ans) <- dim(ans)[1:2]
    ans <- as.data.frame(ans)
    colnames(ans) <- dimnames(x[['results']])[['stat']]
    ans <- data.frame(vars=dimnames(x[['results']])[['vars']], ans)
    row.names(ans) <- 1L:nrow(ans)
    xbalvars <- c("vars", "Control", "Treatment",
                  "adj.diff", "std.diff", "pooled.sd",
                  "z", "p")
    names(xbalvars) <- xbalvars
    vars <- c(xbalvars[setdiff(names(xbalvars), names(varnames_crosswalk))],
              varnames_crosswalk)
    reportme <- intersect(names(vars),
                      c('vars', dimnames(x[['results']])[['stat']])
                      )#to adjust order of columns
    ans <- ans[reportme]
    colnames(ans) <- vars[colnames(ans)]
    ans
}
#' @return
#' @rdname tidy.xbal
#' @S3method
glance.xbal <- function(x, strata=dimnames(x[['results']])[['strata']][1])
{
    ans <- x[['overall']][strata,,drop=FALSE]
    ans <- as.data.frame(ans)
    ans
    }
