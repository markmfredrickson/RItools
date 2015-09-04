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

    theresults <- x$results

    # for historical reasons, what the user requests and the column names in the per-variable table are not the same
    lookup <- c("std.diffs" = "std.diff", "z.scores" = "z", 
                "adj.mean.diffs" = "adj.diff",
                "adj.mean.diffs.null.sd" = "adj.diff.null.sd", "p.values" = "p")

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
      ftable(data, col.vars=c("strata","stat"),row.vars=c("vars"))
    }

    thevartab <- ftabler(theresults) # we'll update this variable later if we include p-values or significance stars

    if (show.signif.stars && !show.pvals && !is.null(theresults) && hasP ) {

      Signif <- signifier(theresults[,"p",which.strata,drop=FALSE])

      ##Nicer alignment, but not as nice labels
      ##junk<-do.call(cbind,lapply(which.strata,function(x){cbind(as.data.frame(theresults[,,x])," "=format(Signif[,,x]))}))
      newresults <- array(dim=dim(theresults)+c(0,1,0),
                          dimnames=list(vars=dimnames(x$results)[["vars"]],
                              stat=c(dimnames(x$results)[["stat"]],"sig."),
                              strata=dimnames(x$results)[["strata"]]))

      newresults[,-grep("sig.", dimnames(newresults)[[2]]),] <- format(theresults,DIGITS)
      newresults[dimnames(Signif)[["vars"]], "sig.",dimnames(Signif)[["strata"]]]<-format(Signif)


      if (horizontal){
        theftab <- ftabler(
            newresults[which.vars, c(which.stats[!(which.stats=="p")],"sig."), which.strata, drop=FALSE])

        attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="sig."] <- ""
        if ("z" %in% which.stats) {
          attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="z"]<-"   z   "
        }

        thevartab<-theftab
      } else {
        thevartab <- sapply(which.strata,
                            simplify=FALSE,
                            function(x) {
                              cbind(
                                  as.data.frame(theresults[which.vars,c(which.stats[!(which.stats=="p")]),x]),
                                  " " = format(Signif[which.vars,,x]))
                            })
      }
    }

    if (show.pvals && hasP && !is.null(theresults)) {
      if (horizontal) {
        theftab <- ftabler(
            theresults[which.vars,which.stats,which.strata,drop=FALSE])
        thevartab <- theftab
      } else {
        thevartab <- sapply(
            which.strata,
            simplify=FALSE,
            function(x) {
              as.data.frame(theresults[which.vars,which.stats,x])
            })
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
        ChiSignif <- signifier(theoverall[which.strata,"p.value"])

        theoveralltab <- cbind(format(theoverall[which.strata,],digits=DIGITS),format(ChiSignif))
        names(theoveralltab)[4]<-" "
      }
      theoveralltab<-format(theoverall[which.strata,],digits=DIGITS)
    } else {
      theoveralltab<-NULL
    }

    if (printme) {
      ## RItools:::print.ftable(thevartab,justify.labels="center",justify.data="right") ##doesn't seem to help the alignment problem
      if (!is.null(theresults)){ print(thevartab) }
      if (!is.null(theoverall) && print.overall) {
        cat("---Overall Test---\n")
        print(theoveralltab)
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
