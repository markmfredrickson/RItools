##' A \code{print} method for balance test objects produced by \code{xBalance} and \code{balanceTest}.
##'
##' @title Printing xBalance and balanceTest Objects
##' @param x An object of class "xbal" which is the result of a call
##'   to \code{xBalance}.
##' @param which.strata The stratification candidates to include in
##'   the printout. Default is all.
##' @param which.stats The test statistics to include. Default is all
##'   those requested from the call to \code{xBalance}.
##' @param which.vars The variables for which test information should
##'   be displayed. Default is all.
##' @param print.overall Should the omnibus test be reported? Default
##'   is \code{TRUE}.
##' @param digits To how many digits should the results be displayed?
##'   Default is \code{max(2,getOptions("digits")-4)}.
##' @param printme Print the table to the console? Default is
##'   \code{TRUE}.
##' @param show.signif.stars Use stars to indicate z-statistics larger
##'   than conventional thresholds. Default is \code{TRUE}.
##' @param show.pvals Instead of stars, use p-values to summarize the
##'   information in the z-statistics. Default is \code{FALSE}.
##' @param horizontal Display the results for different candidate
##'   stratifications side-by-side (Default, \code{TRUE}), or as a
##'   list for each stratification (\code{FALSE}).
##' @param report What to report.
##' @param ... Other arguements. Not currently used.
##' @return \describe{
##' \item{vartable}{The formatted table of variable-by-variable
##'   statistics for each stratification.}
##' \item{overalltable}{If the overall Chi-squared statistic is
##'   requested, a formatted version of that table is returned.}
##' }
##' @seealso \code{\link{xBalance}}, \code{\link{balanceTest}}
##' @export
##' @aliases print print.balancetest
##' @keywords print
##' @examples
##' data(nuclearplants)
##'
##'
##' xb1 <- balanceTest(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data = nuclearplants)
##'
##' print(xb1)
##'
##' print(xb1, show.pvals = TRUE)
##'
##' print(xb1, horizontal = FALSE)
##'
##' ## The following doesn't work yet.
##' \dontrun{print(xb1, which.vars=c("date","t1"),
##'          which.stats=c("adj.means","z.scores","p.values"))}
##'
##' ## The following example prints the adjusted means
##' ## labeled as "treatmentvar=0" and "treatmentvar=1" using the
##' ## formula provided to xBalance().
##'
##' # This is erroring with the change to devtools, FIXME
##' \dontrun{print(xb1,
##'       which.vars = c("date", "t1"),
##'       which.stats = c("pr=0", "pr=1", "z", "p"))}
##'
##' ## Only printing out a specific stratification factor
##' xb2 <- balanceTest(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data = nuclearplants)
##'
##' print(xb2, which.strata = "pt")
print.xbal <- function (x, which.strata=dimnames(x$results)[["strata"]],
                        which.stats=dimnames(x$results)[["stat"]],
                        which.vars=dimnames(x$results)[["vars"]],
                        print.overall=TRUE,
                        digits = NULL, printme=TRUE,
                        show.signif.stars=getOption("show.signif.stars"),
                        show.pvals=!show.signif.stars,
                        horizontal=TRUE,
                        report = NULL, ...) {
  ##Notes: right now we've decided that you can't print both signif stars and p-values. make a choice.

    ## ToDo: check x's include.NA.flags attribute to see if user doesn't want to see them.
    ## If not, apply subset.xbal to get rid of the NotMissing variables.

  # withOptions will allow us to safely reset the digits
  # even if an error is thrown the option should be the same after this function
  DIGITS = ifelse(is.null(digits), max(2, getOption("digits")-4), digits)

  # we'll call this function within a wrapped options block below.
  f <- function() {
    ##makeSigStarsStdNormal <- function(zs) {
    ##  if (length(zs)){##c('','.  ','*  ','** ','***'
    ##    factor(c('','.','*','**','***') [1+
    ##             apply(abs(zs)>=matrix(qnorm(c(.95, .975, .995,.9995)),
    ##                        length(zs),4, byrow=TRUE),1, sum)]
    ##           )} else {character(0)}
    ##                                     }
    if (is.null(report)) {
      report <- attr(x, "report")
    }
    x <- subset(x, vars = which.vars, strata = which.strata, stats = which.stats)
    results_array <- x$results

    # for historical reasons, what the user requests and the column names in the per-variable table are not the same
    lookup <- c("std.diffs" = "std.diff", "z.scores" = "z",
                "adj.mean.diffs" = "adj.diff",
                "p.values" = "p",
                "pooled.sd" = "pooled.sd",
                "adj.mean.diffs.null.sd" = "pooled.sd")

    stopifnot(all(report %in% c("all", "chisquare.test", "adj.means", names(lookup))))

    # adj.means gets expanded into the treated and control group mean columns
    if ("adj.means" %in% report) {
      idx <- report == "adj.means"
      report <- c("Treatment", "Control", report[!idx])
      lookup <- c(Treatment = "Treatment", Control = "Control", lookup)
    }
       
    if (!("all" %in% report)) {
      # on this next line, we use anything in report tha tis also in the names of the lookup table
      # it's a little strange looking, but it does the right thing
      tmp <- lookup[report[report %in% names(lookup)]]

      # likewise, don't grab any columns that aren't there
      results_array <- results_array[, tmp[tmp %in% dimnames(results_array)[["stat"]]], , drop = FALSE]
    }

    ## Mark the columns that will require by-row sigfig handling
    orig_units_columns <- intersect(c("Treatment", "Control", "adj.diff", "pooled.sd"),
                                    dimnames(results_array)[["stat"]])

    hasP <- "p" %in% dimnames(results_array)[["stat"]]

    if("chisquare.test" %in% report || "all" %in% report) { ##Extract the omnibus chisquared test from the xbal object...
      theoverall <- x$overall
    } else {
      theoverall<-NULL ##..or set it to NULL if it does not exist.
    }

    if (length(results_array) == 0 && is.null(theoverall)) {
      stop("There is a problem. Probably all of the variables (",
           all.vars(formula(x)),
           ") are constants within strata. Or else there is some other problem, try debug(RItools:::xBalance) to see what might be going on.")
    }

    if (length(results_array)==0 && !is.null(theoverall)){##The user has requested only the omnibus test and not the tests for the individual variables
      results_array<-NULL
      thevartab<-NULL
    }

    signifier <- function(data) {
      symnum(data,
             corr = FALSE,
             na = FALSE,
             abbr.colnames=FALSE, ##from print.summary.lm
             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
             symbols = c("***", "** ", "*  ", ".  ", "   "))
    }

    ftabler <- function(data) {   ##Summarize the variable-by-variable output array as a flat contingency table
        names(dimnames(data))[names(dimnames(data))=="strata"] <- "strata():"
        ftable(data, col.vars=c("strata","stat"),row.vars=c("vars"))
    }

      if (!is.null(results_array)) {
          results_array_char <- array(dim = dim(results_array), dimnames=dimnames(results_array))
          ## next apply rounding to DIGITS sigfigs, by statistic. If there are
          ## multiple stratifications, then the significant figure position should be set
          ## consistently across them.
          for (rcol in setdiff(dimnames(results_array)[["stat"]], c("p", orig_units_columns)))
          {    #i.e. std.diffs and z stats
              res <- results_array[,rcol,]
              dim(res) <- NULL
              res <- round(res,digits=(DIGITS-1)) 
              results_array_char[,rcol,] <- format(res,digits=DIGITS)
          }
          if ("p" %in% dimnames(results_array)[["stat"]])
          {
              res <- results_array[,"p",]
              dim(res) <- NULL
              results_array_char[,"p",] <- format(res, digits=DIGITS)
          }
          if (!is.null(orig_units_columns))
              results_array_char[,orig_units_columns,] <-
                  original_units_var_formatter(results_array[,orig_units_columns,,drop=FALSE], digits=DIGITS)

          thevartab <- ftabler(results_array_char) # we'll update this variable later if we include p-values or significance stars
    }

      if (show.signif.stars && !show.pvals && !is.null(results_array) && hasP )
      {
        Signif <- signifier(results_array[,"p",,drop=FALSE])
        dimnames(Signif)[['stat']] <- 'sig.'
        results_array_char <- abind(results_array_char, format(Signif),
                                    along=2, use.first.dimnames=TRUE)
        names(dimnames(results_array_char)) <- names(dimnames(results_array))

      if (horizontal){
        tmp <- dimnames(results_array_char)[[2]]
        theftab <- ftabler(results_array_char[, c(tmp[!(tmp == "p")]), , drop=FALSE])

        attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="sig."] <- ""
        if ("z" %in% dimnames(results_array_char)$stat) {
          attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="z"]<-"   z   "
        }

        thevartab<-theftab
      } else {
        tmp <- dimnames(results_array_char)[[2]]
        thevartab <- sapply(dimnames(results_array_char)[[3]],
                            function(s) {
                                tmpdata <- results_array_char[, c(tmp[!(tmp == "p")]), s, drop = FALSE]
                                tmpdn <- dimnames(tmpdata)[1:2]
                                dim(tmpdata) <- dim(tmpdata)[1:2]
                                dimnames(tmpdata) <- tmpdn
                                tmpdata <- as.data.frame(tmpdata)
                              cbind(tmpdata," " = format(Signif[,,s]))
                            },
                            simplify=FALSE, USE.NAMES=TRUE
                            )
      }
    }

    if (show.pvals && hasP && !is.null(results_array)) {
      if (horizontal) {
        theftab <- ftabler(results_array_char)
        thevartab <- theftab
      } else {
        thevartab <- sapply(
            dimnames(results_array_char)[[3]],
            function(x) {
              as.data.frame(results_array_char[,,x])
            },
            simplify=FALSE, USE.NAMES=TRUE
            )
      }
    }

    ##if(show.pvals&!("p"%in% dimnames(results_array)[["stat"]])& !is.null(results_array)) {
    ##  stop("You need to request p-values when calling xBalance.")
    ##} ##irrelevant now.

    if (!is.null(theoverall)) {
      nc <- length(results_array)/2
      latex.annotation <- NULL
      ## paste("\\\\ \\hline Overall",
      ##       paste("\\multicolumn{",nc,"}{c}{",preSig,"}"),
      ##       paste("\\multicolumn{",nc,"}{c}{",postSig,"}"),
      ##       sep=" & ")
      if (show.signif.stars) {
        ChiSignif <- signifier(theoverall[,"p.value"])

        theoveralltab <- cbind(format(theoverall,digits=DIGITS),format(ChiSignif))
        names(theoveralltab)[4]<-" "
      }
      theoveralltab<-format(theoverall,digits=DIGITS)
    } else {
      theoveralltab<-NULL
    }

    if (printme) {
      ## RItools:::print.ftable(thevartab,justify.labels="center",justify.data="right") ##doesn't seem to help the alignment problem
      if (!is.null(results_array)){ print(thevartab) }
      if (!is.null(theoverall) && print.overall) {
        cat("---Overall Test---\n")
        print(theoveralltab, quote=FALSE)
        if (show.signif.stars && !show.pvals && hasP) {
          if (!is.null(results_array)) {
            thelegend<-attr(Signif, "legend") ##if we are showing thevartab use the legend from that object
          }
          if (is.null(results_array) && !is.null(theoverall)) {
            thelegend<-attr(ChiSignif,"legend") ##use legend from the overall object if only showing that one
          }
          cat("---\nSignif. codes: ", thelegend, "\n")
        }
      }
      ##cat(paste("n = ", n, ", k = ", k,
      ##          "\nresidual sd = ", fround(summ$sigma, digits), ", R-Squared = ", fround(summ$r.squared, 2),
      ##          "\n", sep = ""))
    } else {
      list(vartable=thevartab,overalltable=theoveralltab)}
  }
  withOptions(list(digits = DIGITS), f)
}


##' formats a var-stat-strata array by var, with rounding
##' potentially rounding a bit less for an "adj.diff" column
##'
##' 
##' @title Formatting suitable for stat expressed in units specific to var
##' @param arr numeric array
##' @param digits number of digits for rounding
##' @param var_format A list of lists. Each named item of the outer list will be matched to a variable. The inner lists should have two items, `mean` and `diff`. The first formats statistics based on averages. The `diff` item should format statistics that are differences.
##' @return array of same dimension as arr but of type character
##' @author Hansen
##' @keywords internal
original_units_var_formatter <- function(arr, digits, var_format = list())
{
    stopifnot(length(dim(arr))==3, names(dimnames(arr))[1]=="vars" )
    newarr <- array(dim=dim(arr), dimnames=dimnames(arr))
    for (vv in dimnames(arr)[["vars"]])
    {
        res <- arr[vv, ,]
        dim(res) <- NULL
        newarr[vv, ,] <- format(res, digits=digits)
    }

    ## if there's an adj.diff column in addition to a Treatment and a Control column,
    ## permit a little more rounding for the latter than the former.
    if (any(dimnames(arr)[[2]]=="adj.diff") & dim(arr)[2]>1) {
        for (vv in dimnames(arr)[["vars"]]) {
            ouc1 <- setdiff(dimnames(arr)[[2]], "adj.diff")
            res <- arr[vv, ouc1,]
            dim(res) <- NULL
            newarr[vv, ouc1,] <- format(res, digits=digits)
        }
    }


    if (length(var_format) > 0) {
        for (i in seq_along(var_format)) {
            whichvar <- names(var_format)[i]
            if (any(whichvar %in% dimnames(newarr)[["vars"]])) {
                fns <- var_format[[i]]

                available <- intersect(c("Treatment", "Control"), dimnames(newarr)[["stat"]])
                for (j in available) {
                    newarr[whichvar, j, ] <- format(fns$mean(arr[whichvar, j, ]), digits = digits)
                }

                available <- intersect(c("adj.diff", "sd.diff", "pooled.sd"), dimnames(newarr)[["stat"]])
                for (j in available) {
                    newarr[whichvar, j, ] <- format(fns$diff(arr[whichvar, j, ]), digits = digits)
                }
            }
        }
    }

    return(newarr)
}
