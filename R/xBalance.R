xBalance <- function(fmla,
                     data,
                     strata = NULL,
                     report=c("std.diffs","z.scores","adj.means","adj.mean.diffs","adj.mean.diffs.null.sd",
                         "chisquare.test","p.values", "all")[1:2],
                     #                     include.means=FALSE, chisquare.test=FALSE,
                     stratum.weights = harmonic,
                     na.rm = FALSE,
                     impfn = median,
                     covariate.scaling = NULL,
                     normalize.weights = TRUE,
                     post.alignment.transform = NULL) {

  if (!is.null(strata)) {
    stop("The strata argument has been deprecated. Use 'z ~ x1 + x2 + strata(s)' instead. See ?xBalance for details.")
  }

  stopifnot(is.null(post.alignment.transform) || is.function(post.alignment.transform))

  # Using charmatch instead of pmatch to distinguish between no match and ambiguous match. It reports
  # -1 for no match, and 0 for ambiguous (multiple) matches.
  valid.for.report <- c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","chisquare.test",
                                     "std.diffs","z.scores","p.values","all")
  report.good <- charmatch(report, valid.for.report, -1)

  if (any(report.good == -1)) {
    stop(paste("Invalid option(s) for report:", paste(report[report.good == -1], collapse=", ")))
  }
  if (any(report.good == 0)) {
    stop(paste("Option(s) for report match multiple possible values:", paste(report[report.good == 0], collapse=", ")))
  }

  # Now that we've found the partial matches, get their proper names
  report <- valid.for.report[report.good]

  if("all" %in% report)
    report <- c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","chisquare.test", "std.diffs","z.scores","p.values")

  design <- makeDesign(fmla, data)
  design.weighted <- weightedDesign(design, stratum.weights, normalize.weights)

  descriptives <- weightedDesignToDescriptives(design.weighted, covariate.scaling)

  aggDesign             <- aggregateDesign(design)
  aggDesign.weighted    <- weightedDesign(aggDesign, stratum.weights, normalize.weights)
  aggDesign.transformed <- designAlignStrata(aggDesign.weighted, post.align.transform)

  ### Call xBalanceEngine here.

  RES <- lapply(names(ss.df),
                function(nm) {
                  xBalanceEngine(factor(ss.df[,nm]),
                                 zz,
                                 mm1,
                                 report, swt.ls[[nm]],
                                 s.p, normalize.weights,
                                 post.alignment.transform)
                })
  names(RES) <- names(ss.df)
  ##nms <- paste(rep(names(ss.df), rep(length(RES[[1]]$dfr),length(ss.df))),
  ##            names(RES[[1]]$dfr), sep=".")
  ans <- list() ##the overall function still returns a list because of the overall test info.
  ##results is an array of variables by balance statistics by stratification.
  ##here assuming that the variables and statistics are the same across stratifications (including unstratified).
  ans$results<-array(dim=c(vars=nrow(RES[[1]][["dfr"]]),stat=ncol(RES[[1]][["dfr"]]),strata=length(RES)),
                     dimnames=list(vars=rownames(RES[[1]][["dfr"]]),stat=colnames(RES[[1]][["dfr"]]),strata=names(RES)))

  attr(ans$results, "originals") <- design@OriginalVariables

  for (i in names(RES)) {
    ##print(i);print(RES[[i]][["dfr"]])
    ans$results[,,i]<-as.matrix(RES[[i]][["dfr"]])
  }
  ##dimnames(ans)[["stat"]][grep("Tx",dimnames(ans)[["stat"]])]<-c("adj.mean.strata=0","adj.mean.strata=1")
  ##ans$by.variable <- do.call(cbind, lapply(RES, function(x) x[['dfr']]) )
  ##colnames(ans$by.variable) <- nms
  attr(ans, "fmla") <- formula(fmla)

  if ("chisquare.test" %in% report) {
    ans$overall <- data.frame(chisquare = numeric(length(RES)),
                              df        = numeric(length(RES)),
                              p.value   = numeric(length(RES)),
                              row.names = names(RES))
    for (nn in names(RES)) {
      ans$overall[nn,'chisquare'] <- RES[[nn]]$chisq['chisquare']
      ans$overall[nn,'df']        <- RES[[nn]]$chisq['df']
      ans$overall[nn,'p.value']   <- pchisq(RES[[nn]]$chisq['chisquare'],
                                            df = RES[[nn]]$chisq['df'],
                                            lower.tail = FALSE)

    }

    attr(ans$overall, "tcov") <- lapply(RES, function(r) {
      r$tcov
    })
  }
  class(ans) <- c("xbal", "list")
  ans
}

xBalance.make.stratum.mean.matrix <- function(ss, mm) {
  stopifnot(inherits(ss, "factor")) # just in case a numeric variable is passed in.
  
  post.nobs <- dim(mm)[1]
  nlev <- nlevels(ss)

  # for this matrix, finding the indices of the rows is easy, as there is only one
  # item per row, and there post.nobs number of rows.
  tR <- new("matrix.csr",
            ja = as.integer(as.integer(ss)),
            ia = as.integer(1:(post.nobs+ 1)),
            ra = unsplit(1/tapply(ss,ss,length),ss),
            dimension = c(post.nobs,nlev))

  # With many items per row, we need to break the list of strata
  # down into indices of where each row starts and ends
  # e.g. Say ss = 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 the row indices would be
  # 1 5 12 16 (where 16 is the start of the non existant 4th row)

  L <- new("matrix.csr", #ifelse(oldver,"tripletMatrix","dgTMatrix"),
           ia = as.integer(1:(post.nobs + 1)),
           ja = as.integer(as.integer(ss)),
           ra = rep(1,length(ss)),
           dimension = c(post.nobs,nlev))


  msmn <- t(tR) %*% mm
  msmn <- L %*% msmn

  msmn <- as.matrix(msmn)

  return(msmn)
}


# Extract `strata(...)` arguments from a formula.
findStrata <- function(x, data) {

  t <- terms(x, specials = "strata", data = data)

  strata <- rownames(attr(t, "factors"))[attr(t, "specials")$strata]
  if (length(strata) > 0) {
    # Trying to update(x) directly was causing errors about having a "."
    # and no data. Updating the terms returns a fmla and bypasses the bug.
    x <- update(terms(x, data=data),
                as.formula(paste("~ . - ", paste(strata, collapse="-"))))

    # The gsubs return only the `...` inside `strata(...)`
    return(list(newx = x,
                strata = gsub("\\)", "", gsub("strata\\(", "", strata))))
  }

  return(list(newx = x, strata = NULL))
}
