#' `broom::tidy()`/`glance()` methods for `balanceTest()` results
#'
#' Portion out the value of a [RItools::balanceTest()] call in a manner consistent
#' with assumptions of the broom package.  
#' 
#' [RItools::tidy.xbal()] gives per-variable
#' statistics whereas [RItools::glance.xbal()] extracts combined-difference related
#' calculations. In both cases one has to specify which stratification one wants
#' statistics about, as xbal objects can store info about several stratifications.
#' [RItools::tidy.xbal()] has a parameter `varnames_crosswalk` not shared with
#' [RItools::glance.xbal()]. It should be a named character vector, the elements
#' of which give names of columns to be returned and the names of which correspond
#' to columns of xbal objects' \sQuote{results} entry.  Its ordering dictates the order
#' of the result. The default value translates between conventional xbal
#' column names and broom package conventional names.
#'
#' \describe{
#'     \item{vars}{variable name}
#'     \item{Control}{mean of LHS variable = 0 group}
#'     \item{Treatment}{ mean of LHS variable = 1 group}
#'     \item{adj.diff}{T - C diff w/ direct standardization for strata if applicable}
#'     \item{std.diff}{adj.diff/pooled.sd}
#'     \item{pooled.sd}{pooled SD}
#'     \item{statistic}{`z` column from the xbal object}
#'     \item{p.value}{`p` column from the xbal object}
#' }
#' Additional parameters beyond those listed here are ignored (at this time).
#'
#' @param x object of class `"xbal"`, result of [RItools::balanceTest()]
#'    or [RItools::xBalance()]
#' @param strata which stratification to return info about? Defaults
#'    to last one specified in originating function call (which appears first in the xbal array).
#' @param varnames_crosswalk character vector of new names for xbal columns, named by the xbal column
#' @param format if true, apply `[RItools:::original_units_var_formatter()]` to suitable sub-array en route
#' @param digits passed to `[RItools:::original_units_var_formatter()]`
#' @param ... Additional arguments passed to `[RItools:::original_units_var_formatter()]`
#' @return data frame composed of: for `[RItools::tidy()]`, a column of variable labels (`vars`) and 
#'         additional columns of balance-related stats; for `[RItools::glance()]`, scalars describing 
#'         a combined differences test, if found, and otherwise `NULL`.
#' @export
#' @md
tidy.xbal <- function(x,
                      strata = dimnames(x[['results']])[['strata']][1],
                      varnames_crosswalk = c("z"="statistic", "p"="p.value"),
                      format = FALSE,
                      digits=max(2, getOption("digits")-4),...
                      )
{
    ans <- x[['results']][,,strata, drop=FALSE]
    dim(ans) <- dim(ans)[1:2]
    ans <- as.data.frame(ans)
    colnames(ans) <- dimnames(x[['results']])[['stat']]

    orig_units_columns <- intersect(c("Treatment", "Control", "adj.diff", "pooled.sd"),
                                    dimnames(x[['results']])[['stat']])
    if (format && length(orig_units_columns))
    {
        fres <- original_units_var_formatter(x[['results']][,orig_units_columns,,drop=FALSE], digits=digits, ...)
        fres <- fres[,,strata, drop=FALSE]
        dim(fres) <- dim(fres)[1:2]
        ans[,orig_units_columns] <- fres
        }
    
    ans <- data.frame(vars=dimnames(x[['results']])[['vars']], ans)
    ans$NA.info <- if (is.null(attr(x[['results']], "NMpatterns")))
                   {
                       NA_real_
                   } else attr(x[['results']], "NMpatterns")
    row.names(ans) <- 1L:nrow(ans)
    xbalvars <- c("vars", "Control", "Treatment",
                  "adj.diff", "std.diff", "pooled.sd",
                  "z", "p", "NA.info")
    names(xbalvars) <- xbalvars
    vars <- c(xbalvars[setdiff(names(xbalvars), names(varnames_crosswalk))],
              varnames_crosswalk)
    reportme <- colnames(ans) %in% names(vars)
    ans <- ans[reportme]
    colnames(ans) <- vars[colnames(ans)]
    ans
}
#' @rdname tidy.xbal
#' @export
glance.xbal <- function(x, strata=dimnames(x[['results']])[['strata']][1], ...)
{
    ans <- x[['overall']][strata,,drop=FALSE]
    ans <- as.data.frame(ans)
    ans
    }
