xBalance <- function(fmla,
                     data,
                     strata = NULL,
                     report=c("std.diffs","z.scores","adj.means","adj.mean.diffs","adj.mean.diffs.null.sd",
                         "chisquare.test","p.values", "all")[1:2],
                     #                     include.means=FALSE, chisquare.test=FALSE,
                     stratum.weights = harmonic,
                     na.rm = FALSE,
                     impfn = median,
                     covariate.scaling = TRUE,
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

  strataAligned <- alignDesignByStrata(aggDesign.weighted, post.alignment.transform)

  tmp <- lapply(strataAligned, function(i) { do.call(alignedToInferentials, c(list(aggDesign.weighted@Z), i)) })
  names(tmp) <- names(aggDesign@Strata)

  ans <- list()

  # append the z and p to the "descriptives" array (making it somewhat misnamed)
  tmp.z <- as.data.frame(lapply(tmp, function(tt) { tt$z }))
  tmp.p <- as.data.frame(lapply(tmp, function(tt) { tt$p }))
  nstats.previous <- dim(descriptives)[2]
  descriptives <- abind(descriptives, along = 2, tmp.z, tmp.p)
  names(dimnames(descriptives)) <- c("vars", "stat", "strata")

  dimnames(descriptives)[[2]][nstats.previous + 1:2] <- c("z", "p")
  
  inferentials <- do.call(rbind, lapply(tmp, function(s) {
    c(s$csq, s$DF, pchisq(s$csq, df = s$DF, lower.tail = FALSE))
  }))
  colnames(inferentials) <- c("chi.squared", "df", "p.value")
 
  # the meat of our xbal object
  ans$overall <- inferentials
  ans$results <- descriptives

  attr(ans$results, "originals") <- design@OriginalVariables
  attr(ans$overall, "tcov") <- lapply(tmp, function(r) {
    r$tcov
  })
  attr(ans, "fmla") <- formula(fmla)
  attr(ans, "report") <- report # hinting to our summary method later

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
