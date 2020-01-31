##' A \code{print} method for \code{xBalance} objects.
##'
##' @title Printing xBalance Objects
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
##' @seealso xBalance
##' @export
##' @aliases print
##' @keywords print
##' @examples
##' data(nuclearplants)
##'
##' xb0 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'               data=nuclearplants)
##'
##' print(xb0)
##'
##' xb1 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          strata = list(unstratified = NULL, pt = ~pt),
##'          data = nuclearplants,
##'          report = c("all"))
##'
##' str(xb1)
##'
##' print(xb1)
##'
##' print(xb1, show.pvals = TRUE)
##'
##' print(xb1, horizontal = FALSE)
##'
##' ## Now, not asking for the omnibus test
##'
##' print(xb1, which.strata = "pt")
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
    theresults <- x$results

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
      theresults <- theresults[, tmp[tmp %in% dimnames(theresults)[["stat"]]], , drop = FALSE]
    }

    ## Mark the columns that will require by-row sigfig handling
    orig_units_columns <- intersect(c("Treatment", "Control", "adj.diff", "pooled.sd"),
                                    dimnames(theresults)[["stat"]])

    hasP <- "p" %in% dimnames(theresults)[["stat"]]

    if("chisquare.test" %in% report || "all" %in% report) { ##Extract the omnibus chisquared test from the xbal object...
      theoverall <- x$overall
    } else {
      theoverall<-NULL ##..or set it to NULL if it does not exist.
    }

    if (length(theresults) == 0 && is.null(theoverall)) {
      stop("There is a problem. Probably all of the variables (",
           all.vars(formula(x)),
           ") are constants within strata. Or else there is some other problem, try debug(RItools:::xBalance) to see what might be going on.")
    }

    if (length(theresults)==0 && !is.null(theoverall)){##The user has requested only the omnibus test and not the tests for the individual variables
      theresults<-NULL
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

    if (!is.null(theresults)) {
    thevartab <- ftabler(theresults) # we'll update this variable later if we include p-values or significance stars
    }

    if (show.signif.stars && !show.pvals && !is.null(theresults) && hasP ) {

      Signif <- signifier(theresults[,"p",,drop=FALSE])

      ##Nicer alignment, but not as nice labels
      ##junk<-do.call(cbind,lapply(which.strata,function(x){cbind(as.data.frame(theresults[,,x])," "=format(Signif[,,x]))}))
      newresults <- array(dim = dim(theresults) + c(0,1,0),
                          dimnames=list(vars=dimnames(theresults)[["vars"]],
                              stat=c(dimnames(theresults)[["stat"]],"sig."),
                              strata=dimnames(theresults)[["strata"]]))
        ## next apply rounding to DIGITS sigfigs, by statistic. If there are
        ## multiple stratifications, then the significant figure position should be set
        ## consistently across them.
        for (rcol in setdiff(dimnames(theresults)[["stat"]], orig_units_columns))
        {
            res <- theresults[,rcol,]
            dim(res) <- NULL
            if (rcol!="p") { #std.diffs and z stats
                res <- round(res,digits=(DIGITS-1)) }
            newresults[,rcol,] <- format(res,digits=DIGITS)
        }
        if (!is.null(orig_units_columns))
            for (vv in dimnames(theresults)[["vars"]])
            {
                res <- theresults[vv, orig_units_columns,]
                dim(res) <- NULL
                newresults[vv, orig_units_columns,] <- format(res, digits=DIGITS)
            }
        ## if there's an adj.diff column in addition to a Treatment and a Control column,
        ## permit a little more rounding for the latter than the former.
        if (any(orig_units_columns=="adj.diff") & length(orig_units_columns)>1)
                        for (vv in dimnames(theresults)[["vars"]])
                        {
                            ouc1 <- setdiff(orig_units_columns, "adj.diff")
                            res <- theresults[vv, ouc1,]
                            dim(res) <- NULL
                            newresults[vv, ouc1,] <- format(res, digits=DIGITS)
                        }

            
      newresults[dimnames(Signif)[["vars"]], "sig.",dimnames(Signif)[["strata"]]]<-format(Signif)


      if (horizontal){
        tmp <- dimnames(newresults)[[2]]
        theftab <- ftabler(newresults[, c(tmp[!(tmp == "p")]), , drop=FALSE])

        attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="sig."] <- ""
        if ("z" %in% dimnames(newresults)$stat) {
          attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="z"]<-"   z   "
        }

        thevartab<-theftab
      } else {
        tmp <- dimnames(newresults)[[2]]
        thevartab <- sapply(dimnames(newresults)[[3]],
                            function(s) {
                                tmpdata <- newresults[, c(tmp[!(tmp == "p")]), s, drop = FALSE]
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

    if (show.pvals && hasP && !is.null(theresults)) {
      if (horizontal) {
        theftab <- ftabler(theresults)
        thevartab <- theftab
      } else {
        thevartab <- sapply(
            dimnames(theresults)[[3]],
            function(x) {
              as.data.frame(theresults[,,x])
            },
            simplify=FALSE, USE.NAMES=TRUE
            )
      }
    }

    ##if(show.pvals&!("p"%in% dimnames(theresults)[["stat"]])& !is.null(theresults)) {
    ##  stop("You need to request p-values when calling xBalance.")
    ##} ##irrelevant now.

    if (!is.null(theoverall)) {
      nc <- length(theresults)/2
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
      if (!is.null(theresults)){ print(thevartab) }
      if (!is.null(theoverall) && print.overall) {
        cat("---Overall Test---\n")
        print(theoveralltab, quote=FALSE)
        if (show.signif.stars && !show.pvals && hasP) {
          if (!is.null(theresults)) {
            thelegend<-attr(Signif, "legend") ##if we are showing thevartab use the legend from that object
          }
          if (is.null(theresults) && !is.null(theoverall)) {
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
