###############################################################################
# ModelMatrixPlus Objects: covariates with term-specific missingness info & indexing
################################################################################

setClassUnion("Contrasts", c("list", "NULL"))
##' ModelMatrixPlus S4 class
##'
##' If the Covariates matrix has an intercept, it will only be in the first column.
##'
##' More on NotMissing slot: It's matrix of numbers in [0,1].
##' First col is entirely \code{TRUE} or 1, like an intercept, unless corresponding
##' UnitWeight is 0, in which case it may also be 0 (see below). Subsequent cols
##' present only if there are missing covariate values, in which case these cols are
##' named for terms (of the original calling formula or data frame) that possess
##' missing values.  Terms with the same missing data pattern are mapped to a single
##' column of this matrix.  If the ModelMatrixPlus is representing elements, each column should
##' be all 1s and 0s, indicating which elements have non-missing values for the term
##' represented by that column.  If the ModelMatrixPlus as a whole represents clusters,
##' then there can be fractional values, but that situation should only arise in the
##' DesignOptions class exension of this class, so it's documented there.
##'
##' @slot Covariates The numeric matrix that `model.matrix` would have returned.
##' @slot OriginalVariables look-up table associating Covariates columns with terms of the originating model formula
##' @slot TermLabels labels of terms of the originating model formula
##' @slot contrasts Contrasts, a list of contrasts or NULL, as returned by `model.matrix.default`
##' @slot NotMissing Matrix of numbers in [0,1] with as many rows as the Covariates table but only one more col than there are distinct covariate missingness patterns (at least 1, nothing missing). First col is entirely T or 1, like an intercept.
##' @slot NM.Covariates integer look-up table mapping Covariates columns to columns of NotMissing.  (If nothing missing for that column, this is 0.)
##' @slot NM.terms integer look-up table mapping term labels to columns of NotMissing (0 means nothing missing in that column)
##' @slot UnitWeights vector of weights associated w/ rows of the ModelMatrixPlus
##' @keywords internal
setClass("ModelMatrixPlus",
         slots=c(Covariates="matrix",
                  OriginalVariables="integer",
                  TermLabels="character",
                  Contrasts="Contrasts",
                  NotMissing="matrix",
                  NM.Covariates="integer",
                  NM.terms="integer",
                  UnitWeights = "numeric" )
         )

#' @method as.matrix ModelMatrixPlus
#' @export
as.matrix.ModelMatrixPlus <- function(x, ...)
    {
        ans <- x@Covariates
        attr(ans, "assign") <- x@OriginalVariables
###        attr(ans, "term.labels") <- # a model.matrix wouldn't really
###            x@TermLabels # have this, although maybe it should
        attr(ans, "contrasts") <- x@Contrasts
        ans
    }

##' Grow a model matrix while at the same time compactly
##' encoding missingness patterns in RHS variables of a model frame.
##'
##'
##' @title Model matrices along with compact encodings of data availability/missingness
##' @param object Model formula or terms object (as in `model.matrix`)
##' @param data data.frame, as in `model.matrix()` but has to have \sQuote{\code{(weights)}} column
##' @param remove.intercept logical
##' @param ... passed to `model.matrix.default` (and further)
##' @return ModelMatrixPlus, i.e. model matrix enriched with missing data info
##' @author Ben B Hansen
##' @keywords internal
##' 
model_matrix <- function(object, data = environment(object), remove.intercept=TRUE, ...) {
  # mf <- model.frame(object, data, na.action = na.pass)
  tms <- terms(object)
  term.labels <- attr(tms, "term.labels")

  uweights <- as.vector(model.weights(data))
  if (is.null(uweights))
    stop("model_matrix() expects its data arg to be a model frame containing weights")
  stopifnot(is.numeric(uweights), all(!is.na(uweights)), all(uweights>=0))
  
  covariates <- model.matrix(object = object, data = data, ...) 

  assign <- attr(covariates, "assign")
  attr(covariates, "assign") <- NULL
  stopifnot(all(assign %in% 0:length(term.labels)))

  contrasts <- attr(covariates, "contrasts")
  attr(covariates, "contrasts") <- NULL

  if (remove.intercept &&
      (iloc <- match("(Intercept)", colnames(covariates), nomatch=0))
      ) {
    assign <- assign[-iloc]
    covariates <- covariates[,-iloc, drop=FALSE]
  }

    logicalsT <- (paste0(term.labels, "TRUE") %in% colnames(covariates)) &
        !(paste0(term.labels, "TRUE") %in% term.labels)
    logicalsF <- (paste0(term.labels, "FALSE") %in% colnames(covariates)) &
        !(paste0(term.labels, "FALSE") %in% term.labels)
    if (any(logicalsF & logicalsT))
    {
        rlocs <- colnames(covariates) %in% paste0(term.labels[logicalsF & logicalsT], "FALSE")
        assign <- assign[!rlocs]
        covariates <- covariates[,!rlocs, drop=FALSE]
    }
    if (any(logicalsT))
    {
        slocs <- colnames(covariates) %in% paste0(term.labels[logicalsT], "TRUE")
        colnames(covariates)[slocs] <-
            substr(colnames(covariates)[slocs], 1L,
                   nchar(colnames(covariates)[slocs]) - 4L)
        }

  cols.by.term <- lapply(1L:length(term.labels),
                         function(whichterm) (1:ncol(covariates))[assign==whichterm])
  ccs.by.term <-
    lapply(cols.by.term, function(whichcols) {
      complete.cases(covariates[,whichcols])
    })
  names(ccs.by.term) <- term.labels
  null.record <- rowSums(as.data.frame(ccs.by.term))==0
  if (any(null.record)) #in this case make sure it's NA in each covariate column
      null.record[null.record] <- apply(is.na(covariates[null.record,,drop=FALSE]), 1, all)

  terms.with.missings <- !sapply(ccs.by.term, all)
    
  ccs.by.term <- c(list('_any Xs recorded_'=!null.record),
                   ccs.by.term)
  terms.with.missings <- c(TRUE, # in order always to have same leading entry
                           terms.with.missings)

  nm.covs <- integer(ncol(covariates))
  nm.terms <- integer(length(terms.with.missings))

  nmcols <- ccs.by.term[terms.with.missings]
  nmcols.dupes <- duplicated(nmcols, fromLast=FALSE)
  if (any(nmcols.dupes))
  {
      nm.terms[terms.with.missings][!nmcols.dupes] <- 1L:sum(!nmcols.dupes)
      nm.terms[terms.with.missings][nmcols.dupes] <-
        match(nmcols[nmcols.dupes], nmcols[!nmcols.dupes])
      ## At this point nm.terms is a look-up table with an entry for
      ## _any Xs recorded_ and for each provided term. Entries are
      ## 0 if there are no null records, or no missings on the term in question;
      ## otherwise if there are null records and/or term missings, the
      ## entry gives the column of the matrix `notmissing`, formed a few lines
      ## down, that describes the relevant missingness pattern.
  } else nm.terms[terms.with.missings] <- 1L:length(nmcols)

    notmissing <- as.data.frame(nmcols[!nmcols.dupes])
    names(notmissing) <- names(ccs.by.term)[terms.with.missings][!nmcols.dupes]
    ## Now form the look-up table associating columns of the covariates matrix
    ## with columns of matrix `notmissing` describing relevant missingness patterns.
    ## Columns associated with terms on which there was no missingness get a 0 here.
    nm.terms <- nm.terms[-1L]
    nm.covs[assign>0] <- nm.terms[assign[assign>0]]

    notmissing <- as.matrix(notmissing)
    
  new("ModelMatrixPlus",
      Covariates=covariates,
      OriginalVariables=assign,
      TermLabels=term.labels,
      Contrasts=contrasts,
      NotMissing=notmissing,
      NM.Covariates=nm.covs,
      NM.terms=nm.terms,
      UnitWeights = uweights )
}


################################################################################
# DesignOptions Objects: communicate cluster, strata, treatment, and covariate information
################################################################################

#' DesignOptions S4 class
#'
#' Extends the ModelMatrixPlus class
#'
##' If the DesignOptions represents clusters of elements, as when it was created
##' by aggregating another DesignOptions or ModelMatrixPlus object, then its Covariates
##' and NotMissing slots are populated with (weighted) averages, not totals.  E.g.,
##' NotMissing columns consist of weighted averages of element-wise non-missingness indicators
##' over clusters, with weights given by (the element-level precursor to) the UnitWeights
##' vector.  As otherwise, columns of the NotMissing matrix represent terms
##' from a model formula, rather than columns the terms may have expanded to.
##'
##' If present, the null stratification (all units in same stratum) in indicated
##' by the corresponding column of the \code{StrataFrame} slot bearing the name
##' \sQuote{\code{--}}.

#' @slot Z Logical indicating treatment assignment
#' @slot StrataFrame Factors indicating strata
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @keywords internal
#'
setClass("DesignOptions",
         slots = c(
           Z                 = "logical",
           StrataFrame       = "data.frame",  
             Cluster           = "factor"),
         contains = "ModelMatrixPlus"
         )

#
#
##' Create a DesignOptions object from a formula and some data
##'
##' The formula must have a left hand side that can be converted to a logical.
##'
##' On the RHS:
##'  - It may have at most one cluster() argument.
##'  - It may have one or more strata() arguments.
##'  - All other variables are considered covariates.
##'
##' NAs in a cluster() or strata() variable will be dropped.
##' NAs in covariates will be passed through, but without
##' being flagged as NotMissing (as available data items will)
##'
##' @param fmla Formula
##' @param data Data
##' @return DesignOptions
##' @import stats
##' @keywords internal
##'
makeDesigns <- function(fmla, data) {

    ts <- terms(fmla, data = data[setdiff(colnames(data), '(weights)')],
                specials = c("cluster", "strata"))

  if (attr(ts, "response") == 0) {
    stop("You must include a treatment assignment variable on the left hand side")
  }

  includeUnstratified <- attr(ts, "intercept") == 1
  if (is.null(attr(ts, "specials")$strata) && !includeUnstratified) {
    stop("If you remove the unadjusted comparison (using '-1'), you must include a stratification")
  }

  if (length(attr(ts, "specials")$cluster) > 1) {
    stop("At most one cluster variable is allowed")
  }

  ## For the moment: it is not an error to fail to include at least one covariate.
  ## This would probably be an error further down the chain for most of our applications,
  ## for now, we're keeping this function agnostic on the subject

  # split up the formula into a structural formula (treatment, strata, cluster)
  # and a data component.
  vnames <- rownames(attr(ts, "factors"))
  treatment.name <- vnames[attr(ts, "response")]
  str.vnames <- vnames[c(attr(ts, "specials")$cluster, attr(ts, "specials")$strata)]

  str.fmla <- formula(paste0("factor(", treatment.name, ")", " ~ ", paste0(collapse = "+", c(1, str.vnames))),
                      env=environment(fmla))
  str.tms  <- terms(str.fmla, data = data,
                    specials = c("cluster", "strata"))
  str.data <- model.frame(str.tms, data = data, na.action = na.pass, drop.unused.levels=TRUE)


  ## check that strata and clusters have the proper relationships with treatment assignment
  treatmentCol <- colnames(str.data)[attr(str.tms, "response")]
  clusterCol <- colnames(str.data)[attr(str.tms, "specials")$cluster]
  strataCols <- colnames(str.data)[attr(str.tms, "specials")$strata]

  if (includeUnstratified) {
    str.data$`--` <- 1
    strataCols <- c(strataCols, "--")
  }

  ## first we check cluster info
  if (length(clusterCol) > 0) {

    # each cluster should have 1 and only 1 type of unit (treated or control)
    tbl <- table(str.data[, c(treatmentCol, clusterCol)])
    isGood <- apply(tbl, 2, function(x) { sum(x != 0)  == 1})
    if (!(all(isGood))) {
      stop("The following clusters did not have identical treatment assignment: ",
           paste(collapse = ", ", colnames(tbl)[!isGood]))
    }

    # each cluster should be nested entirely in each strata.
    for (s in strataCols) {
      tbl <- table(str.data[, c(s, clusterCol)], useNA = "ifany")
      isGood <- apply(tbl, 2, function(x) { sum(x != 0) == 1 })
      if (!(all(isGood))) {
        stop("In ", s, ", the following clusters were not nested within stratum levels: ",
             paste(collapse = ", ", colnames(tbl)[!isGood]))
      }
    }

  }

  ## now check that strata works (apart from treatment)
  if (nlevels(str.data[, treatmentCol]) != 2) {
    stop("Treatment assignment must have exactly two levels.")
  }

  ## for each strata, we want there to be at least one treated and control unit

  for (s in strataCols) {

    if (all(is.na(str.data[, s]))) {
      stop("All levels in ", s, " are NA")
    }

    tbl <- table(str.data[, c(treatmentCol, s)])
    isGood <- apply(tbl, 2, function(x) { sum(x != 0) == 2 })
    if (!any(isGood))
        stop("In ", s, " there were no strata with both treatment and control units")
    ## otherwise, just ignore the bad strata

    if (!(all(isGood))) {
        warning(sum(!isGood), " levels of ", s, "  did not include both treated and control units")
    }
  }

  ## OK! data looks good. Let's proceed to make the design object with covariate data

  data.fmla <- update(ts, paste("~", paste0(collapse = " - ", c(".", str.vnames))))
  data.data <- model.frame(data.fmla, data, na.action = na.pass) #
  data.data$'(weights)' <- data$'(weights)'
  
  # Convert remaining character columns to factors so make behavior consistent (issue #127)
  chrs <- colnames(data.data)[sapply(data.data, is.character)]
  for (f in chrs) {
    data.data[, f] <- as.factor(data.data[, f])
  }
  
  # knock out any levels that are not used
  fcts <- colnames(data.data)[sapply(data.data, is.factor)]
  for (f in fcts) {
    data.data[, f] <- factor(data.data[, f])
  }

  # we want our own contrast function for the factors that expands each level to its own dummy
  clist <- lapply(data.data, function(x) {
    if (is.factor(x)) {
      structure(diag(nlevels(x)), dimnames = list(levels(x), levels(x)))
    } else {
      NULL
    }
})
    names(clist) <- colnames(data.data)
    clist <- clist[!sapply(clist, is.null)]

    desm <- model_matrix(terms(data.data), data.data, remove.intercept=TRUE, contrasts.arg = clist)

  if (length(clusterCol) > 0) {
    Cluster <- str.data[, clusterCol]
  } else {
    Cluster <- 1:(dim(str.data)[1])
  }

  Z <- str.data[, treatmentCol]
  tmp <- str.data[, strataCols, drop = FALSE]
  colnames(tmp) <- gsub(colnames(tmp), pattern = "strata\\((.*)\\)", replacement = "\\1")
  strata.frame <- data.frame(lapply(tmp, factor), check.names = FALSE)

  return(new("DesignOptions",
             Z                 = as.logical(as.numeric(Z) - 1), #b/c it was built as factor
             StrataFrame       = strata.frame,
             Cluster           = factor(Cluster),
             UnitWeights = desm@UnitWeights,
             Covariates = desm@Covariates,
                  OriginalVariables=desm@OriginalVariables,
                  TermLabels=desm@TermLabels,
                  Contrasts=desm@Contrasts,
                  NotMissing=desm@NotMissing,
                  NM.Covariates=desm@NM.Covariates,
             NM.terms=desm@NM.terms)
         )
}

# Add stratum weights to a design
#' Stratum Weighted DesignOptions
#'
#' @slot Sweights stratum weights
setClass("StratumWeightedDesignOptions",
         slots = list(Sweights = "list"),
         contains = "DesignOptions")


##' Create stratum weights to be associated with a DesignOptions
##'
##' Apply weighting function to a DesignOptions by stratum, returning results in a format
##' suitable to be associated with an existing design and used in further calcs.
##'
##' The function expects its DesignOptions argument to represent aggregated data,
##' i.e. clusters not elements within clusters.  Its \code{stratum.weights} argument
##' is a function that is applied to a data frame representing clusters,
##' with variables \code{Tx.grp}, \code{stratum.code}, covariates as named in the
##' \code{design} argument (a DesignOptions object), and \code{unit.weights} (either as
##' culled or inferred from originating \code{\link{balanceTest}} call or as aggregated up
##' from those unit weights).  Returns a
##' weighting factor to be associated with each stratum, this factor determining the stratum
##' weight by being multiplied by mean of unit weights over clusters in that stratum.
##'
##' Specifically, the function's value is a data frame of two variables,
##' \code{sweights}  and \code{wtratio}, with rows representing strata.
##' The \code{sweights} vector represents internally
##' calculated or user-provided \code{stratum.weights}, one for each
##' stratum, scaled so that their sum is 1; in Hansen & Bowers (2008), these
##' weights are denoted \eqn{w_{b}}. \code{wtratio} is the ratio of
##' \code{sweights} to the product of half the harmonic
##' mean of \eqn{n_{tb}} and \eqn{n_{cb}}, the number of treatment and control
##' clusters in stratum \eqn{b}, with the mean of the weights associated with
##' each of these clusters.  In the notation of Hansen & Bowers
##' (2008), this is \eqn{w_{b}/(h_b \bar{m}_b)}. Despite the name
##' \sQuote{\code{wtratio}}, this ratio's denominator is not a weight
##' in the sense of summing to 1 across strata.  The ratio is expected
##' downstream in \code{HB08} (in internal calculations
##' involving \sQuote{\code{wtr}}).
##'
##'
##' @param design DesignOptions
##' @param stratum.weights Stratum weights function. Will be fed a count data.frame with Tx.grp (indicating the treatment group), stratum.code, all other covariates and unit.weights.
##' @return data frame w/ rows for strata, cols \code{sweights} and \code{wtratio}.
##'
##' @keywords internal

DesignWeights <- function(design, stratum.weights = harmonic_times_mean_weight) {
  stopifnot(inherits(design, "DesignOptions"),
            is.function(stratum.weights),
            !is.null(design@UnitWeights),  # TODO: take these checks out if we do not use unit weights in this function
            all(design@UnitWeights >= 0 ) )

  lapply(design@StrataFrame, function(s) {

    stratifier <- factor(s)

    swts <- do.call(stratum.weights,
                    args = list(data =
                                  data.frame(Tx.grp = design@Z,
                                             stratum.code = stratifier,
                                             design@Covariates,
                                             unit.weights=design@UnitWeights,
                               check.names = FALSE)),
                envir=parent.frame())

    ## checking that the returned swts meet their requirements
    if (all(is.na(swts)))
      stop(paste("All stratum weights NA."))

    if (any(is.na(swts))) {
      swts[is.na(swts)] <- 0
      warning(paste("NAs in stratum weights; to be interpreted as 0s.", sep=""))
    }

    if (any(swts<0))
      stop("stratum weights must be nonnegative")

    if (identical(harmonic_times_mean_weight, stratum.weights)) {
        ## Subtlety re correspondence of codebase with H&B08:
        ## H&B's h_b equals *half* the harmonic mean of n_t, n_c;
        ## RItools's harmonic() calculates harmonic means of
        ## (n_t, n_c) pairs, w/o the (1/2) factor. (Likewise for
        ## harmonic_times_mean_weight().)
        ## If we're here, then the current sweights
        ## will have been calculated as (2 * h_b * m-bar_b), by
        ## harmonic_times_mean_weight(). Since wtratio needs to
        ## compare these weights to (h_b * m-bar_b) as in H&B, we have:
        wtratio <- rep(2, length(swts))
        ## (normalization of sweights to be addressed below).
    } else {
        hwts <- harmonic_times_mean_weight(
            data.frame(Tx.grp = design@Z,
                       stratum.code=stratifier,
                       unit.weights=design@UnitWeights,
                       check.names = FALSE))
      wtratio <- swts/(hwts/2) # Re 1/2 factor, see note immediately above
    }

    sweights <- swts
    sum.sweights <- sum(sweights, na.rm=TRUE)
    sweights <- sweights / sum.sweights

    wtratio <- wtratio / sum.sweights
    names(wtratio) <- levels(stratifier)

    data.frame(sweights = sweights,
               wtratio = wtratio,
               row.names = levels(stratifier))
  })
}

##' Generate Descriptives
##'
##' Use a design object to generate descriptive statistics that ignore clustering.
##' Stratum weights are respected if provided (by passing a design arg of
##' class StratumWeightedDesignOptions). If not provided, stratum weights
##' default to "Effect of Treatment on Treated" weighting.  That is, when
##' combining within-stratum averages (which will themselves have been
##' weighted by unit weights), each stratum receives a weight equal to
##' the product of the stratum sum of unit weights with the fraction
##' of clusters within the stratum that were assigned to the treatment
##' condition.
##'
##' By default, covariates are scaled by their pooled s.d.s, square roots
##' of half of their treatment group variances plus half of their control
##' group variances.  If weights are provided, these are weighted variances.
##' If descriptives are requested for an unstratified setup, i.e. a 
##' stratification named \sQuote{\code{--}}, then covariate s.d.s
##' are calculated against it; otherwise the variances reflect stratification,
##' and are calculated against the first stratification found.  Either way,
##' if descriptives are calculated for multiple stratifications, only one
##' set of covariate s.d.s will have been calculated, and these underlie
##' standard difference calculations for each of the stratifications.
##'
##' If a named numeric \code{covariate.scales} argument is provided, any
##' covariates named in the vector will have their pooled s.d.s taken from
##' it, rather than from the internal calculation. 
##' 
##' @param design A DesignOptions object
##' @param covariate.scales Scale estimates for covariates, a named numeric vector
##' @return Descriptives
##' @keywords internal
##'
designToDescriptives <- function(design, covariate.scales = NULL) {
  stopifnot(inherits(design, "DesignOptions")) # defensive programming
  covars <- ifelse(is.na(design@Covariates), 0, design@Covariates)

  ## Tack NM cols onto covars, but with intercept col listed last
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  covars <- cbind(covars, design@NotMissing[, NMcolperm])
  vars <- c(colnames(design@Covariates),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
  colnames(covars)   <-  vars
  covars.nmcols <- c(1L + design@NM.Covariates, rep(1L, k.NM ) )

  covariate.scales  <-
      if (is.null(covariate.scales) |
          (bad  <- !is.null(covariate.scales) &
               (  !is.numeric(covariate.scales) | is.null(names(covariate.scales)) )
              )
          )
      {
          if (bad)
              warning("'covariate.scales' should be NULL or a named numeric; ignoring")
          setNames(numeric(0), character (0))
      } else {
          common_names  <- intersect(names(covariate.scales), vars)
          covariate.scales[common_names]
      }  

  stratifications <- colnames(design@StrataFrame)

  Uweights <- design@UnitWeights * cbind(1, design@NotMissing)

  ans <- array(NA,
               dim = c(length(vars), 5, length(stratifications)),
               dimnames = list(
                   "vars" = vars,
                   "stat" = c("Control", "Treatment", "std.diff", "adj.diff", "pooled.sd"),
                   "strata" = stratifications))
  if (any(null_strat  <- stratifications=="--"))
      stratifications  <- c("--", stratifications[!null_strat])
      
  for (s in stratifications) {

    S <- SparseMMFromFactor(design@StrataFrame[[s]])
    stratlevs <- levels(design@StrataFrame[[s]])
    if (inherits(design, "StratumWeightedDesignOptions") & length(stratlevs)>1)
        {
            Sweights <- design@Sweights[[s]][['sweights']]
            if (!all(names(Sweights) %in% stratlevs))
                stop(paste("Strata weights/names mismatch, stratification", s) )
            Swts <- rep(0, length(stratlevs))
            names(Swts) <- stratlevs
            Swts[names(Sweights)] <- Sweights
        } else {
          if (length(stratlevs)==1) {Swts <- 1
          ## if the stratifier has just 1 level, we don't
          ## need stratum weights, so just use 1.  If it has more
          ## than 1 level, then we're here because stratum weights
          ## were not passed down, and we'll need to use defaults,
          ## which are more convenient to set only implictly; see further down.
            names(Swts) <- stratlevs} else Swts <- NULL
                }

    Z <- as.numeric(design@Z)
    ZZ <- S * Z
    WW <- S * (1 - Z)


    X.use  <- covars * Uweights[, covars.nmcols]
    X2.use <- covars^2 * Uweights[, covars.nmcols]

    n1 <- t(Uweights) %*% Z
    n1 <- n1[covars.nmcols, ]
    n0 <- t(Uweights) %*% (1 - Z)
    n0 <- n0[covars.nmcols, ]

    ## Now calculate assignment/stratum weights
    cluster.representatives <- !duplicated(design@Cluster)
    txclus.by.strat <- t(ZZ) %*% as.matrix(cluster.representatives)
    txclus.by.strat <- as.matrix(txclus.by.strat)
    ctlclus.by.strat <- t(WW) %*% as.matrix(cluster.representatives)
    ctlclus.by.strat <- as.matrix(ctlclus.by.strat)
    nclus.by.strat <- txclus.by.strat + ctlclus.by.strat
    strat.sum.uweights <- t(S) %*% as.matrix(design@UnitWeights)
    strat.sum.uweights <- as.matrix(strat.sum.uweights)
    if (!is.null(Swts)) ## this means stratum weights were passed
    {
    ## Horwitz Thompson-type assignment weights
    tx.wt <- ifelse(txclus.by.strat, nclus.by.strat/txclus.by.strat, 0)
    ctl.wt <-ifelse(ctlclus.by.strat, nclus.by.strat/ctlclus.by.strat, 0)
    tx.wt <- ifelse(strat.sum.uweights, tx.wt/strat.sum.uweights, 0)  # approp. HT weights to estimate, by stratum,
    ctl.wt <-ifelse(strat.sum.uweights, ctl.wt/strat.sum.uweights, 0) # (total uweighted measurements)/(total uweights)
    ## Next factor in stratum weights, `Swts`.
    tx.wt <- tx.wt * Swts
    ctl.wt <- ctl.wt * Swts
    } else { # in this condition no Swts were passed, so we default to
      ## strat.sum.uweights *(txclus.by.strat / nclus.by.strat)
      tx.wt <- ifelse(txclus.by.strat, 1, 0)
      ctl.wt <-ifelse(ctlclus.by.strat, txclus.by.strat/ctlclus.by.strat, 0)
      }
    ## now expand up tx.wt and ctl.wt to match dimensions of data
    tx.wts <- as.vector(as.matrix(S %*% tx.wt))
    ctl.wts <- as.vector(as.matrix(S %*% ctl.wt))

    # ok, now that preliminaries are out of the way, compute some useful stuff.
    ## ratio estimates of means for a "domain" equal to intersection of
    ## treatment group with units for which which the covariate is non-missing
    wtsum.tx <-  t(Uweights * tx.wts) %*% Z
    wtsum.tx <- wtsum.tx[covars.nmcols, ]
    treated.avg <- t(X.use *tx.wts) %*% Z / wtsum.tx

    ## ratio estimates of means over similarly defined domains within the control group
    wtsum.ctl <- t(Uweights * ctl.wts) %*% (1 - Z)
    wtsum.ctl <- wtsum.ctl[covars.nmcols, ]
    control.avg <- t(X.use * ctl.wts) %*% (1 - Z) / wtsum.ctl

    ## only perform pooled s.d. calculation the first time.
      if (s==stratifications[1L])
          {
    ## the use of pmax is to correct rounding errors that lead to variances of -1.0e-16
    ## which causes sqrt() to throw a warning.
    var.1 <- (t(X2.use *tx.wts) %*% Z - wtsum.tx * treated.avg^2) / wtsum.tx
    var.1 <- var.1 * ifelse(n1>1, n1/(n1 - 1), 0)
    var.0 <- (t(X2.use * ctl.wts) %*% (1 - Z) - wtsum.ctl * control.avg^2) / wtsum.ctl
    var.0 <- var.0 * ifelse(n0>1, n0/(n0 - 1), 0)
    
    var.1[var.1 < 0] <- 0
    var.0[var.0 < 0] <- 0
    
    pooled <- sqrt((var.1 + var.0) / 2)
    ## if covariate scales were provided, they override what
    ## was just calculated
    if (length(covariate.scales)) {
        pooled[match(names(covariate.scales), rownames(pooled)),]  <- covariate.scales
        }
    }

    adjustedDifference    <- treated.avg - control.avg
    standardizedDifference <- adjustedDifference / pooled

    ans[, , s] <- c(control.avg, treated.avg, standardizedDifference, adjustedDifference, pooled)
  }

  return(ans)
}

##' Aggregate DesignOptions
##'
##' Totals up all the covariates, as well as user-provided unit weights.
##' (What it does to NotMissing entries is described in docs for DesignOptions class.)
##'
##' If \code{design@Cluster} has extraneous (non-represented) levels, they will be dropped.
##'
##' @param design DesignOptions
##' @return another DesignOptions representing the clusters
##' @keywords internal
##'
aggregateDesigns <- function(design) {
  clusters <- factor(design@Cluster)
  n.clusters <- nlevels(clusters)

  if (n.clusters == length(clusters)) {
    return(design)
  }

  dupes <- duplicated(clusters)
  Z <- design@Z[!dupes]
  StrataFrame <- design@StrataFrame[!dupes,, drop = FALSE]
  Cluster <- clusters[!dupes]
  names(Z) <- as.character(Cluster)

  C_transp <- t(SparseMMFromFactor(clusters))
  # To align w/ this `C_transp`, everything to be returned 
  # needs to align w/ `levels(Cluster)`, not with 
  # `Cluster` itself. So,
  cperm <- match(levels(Cluster), as.character(Cluster))
  Z <- Z[cperm]
  StrataFrame  <- StrataFrame[cperm, , drop=FALSE]
  Cluster <- Cluster[cperm]
  
  unit.weights <- as.matrix(C_transp %*% as.matrix(design@UnitWeights))

  dim(unit.weights) <- NULL
  names(unit.weights) <- levels(Cluster)

  Uweights.tall <- design@UnitWeights * design@NotMissing
  Covariates <- as.matrix(C_transp %*% ifelse(Uweights.tall[,pmax(1L,design@NM.Covariates), drop=FALSE],
                                          design@Covariates *
                                              Uweights.tall[,pmax(1L,design@NM.Covariates), drop=FALSE],
                                          0)
                          )
  Uweights <- as.matrix(C_transp %*% Uweights.tall)
  Covariates <- ifelse(Uweights[,pmax(1L,design@NM.Covariates), drop=FALSE] > 0,
                       Covariates/Uweights[,pmax(1L,design@NM.Covariates), drop=FALSE],
                       0)
  NotMissing <- ifelse(matrix(unit.weights>0, nrow(Uweights), ncol(Uweights)),
                       Uweights/unit.weights, 0)
  colnames(NotMissing) <- colnames(design@NotMissing)
  
    Covariates <- as.matrix(Covariates)
    colnames(Covariates)   <- colnames(design@Covariates)
    row.names(Covariates) <- levels(Cluster)

  new("DesignOptions",
      Z = Z,
      StrataFrame = StrataFrame,
      Cluster = Cluster,
      UnitWeights = unit.weights,
      NotMissing = NotMissing,
      Covariates = Covariates,
      OriginalVariables = design@OriginalVariables,
      TermLabels=design@TermLabels,
      Contrasts=design@Contrasts,
      NM.Covariates=design@NM.Covariates,
      NM.terms=design@NM.terms)
}

#' CovsAlignedToADesign S4 class
#'
#' A class for representing covariate matrices after alignment within stratum,
#' for a (single) given stratifying factor.  There can also be a clustering variable,
#' assumed to be nested within the stratifying variable.
#'
#' In contrast to DesignOptions, this class represents the combination of a single Covariates
#' table, realized treatment assignment and treatment assignment scheme, not multiple treatment
#' assignment schemes (designs). These Covariates are assumed to reflect regularization, s.t.
#' missings have been patched with a value and then all of the covariate values have been aligned
#' within each stratum.  In lieu of a NotMissing slot there will be Covariates columns, also
#' centered within a stratum, recording
#' non-missingness of the original data.
#'
#' Ordinarily the StrataWeightRatio slot has an entry for each unit, representing ratio of
#' specified stratum weight to the product of \eqn{h_b} (the harmonic mean of \eqn{n_\{tb\}} and
#' \eqn{n_{cb}}, the counts of treatment and control clusters in stratum \eqn{b}) with \eqn{bar-w_b},
#' (the arithmetic mean of aggregated cluster weights within that stratum). It can also
#' be the numeric vector 1, without names, meaning the intended weight ratio is always 1.
#' 
#' @slot Covariates Numeric matrix, as in ModelMatrixPlus, except: will include NM columns; all columns presumed to have been stratum-centered (aligned)
#' @slot UnitWeights vector of weights associated w/ rows of Covariates
#' @slot Z Logical indicating treatment assignment
#' @slot StrataMatrix A sparse matrix with n rows and s columns, with 1 if the unit is in that stratification
#' @slot StrataWeightRatio For each unit, ratio of stratum weight to \eqn{h_b}; but see Details.
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @slot OriginalVariables Look up table associating Covariates cols to terms in the calling formula, as in ModelMatrixPlus
#' @keywords internal
setClass("CovsAlignedToADesign",
         slots =
             c(Covariates="matrix",
             Z                 = "logical",
             StrataMatrix    = "matrix.csr", 
             StrataWeightRatio = "numeric",
             OriginalVariables="integer",
             Cluster = "factor"
             )
         )
# apply this & pass through en route to svd
#' Scale DesignOptions
#' @method scale DesignOptions
##' @param x DesignOptions object
##' @param center logical, or a function acceptable as \code{post.alignment.transform} arg of \code{alignDesignsByStrata()}
##' @param scale logical, whether to scale
scale.DesignOptions  <- function(x, center=TRUE, scale=TRUE)
{
    stopifnot(is(x, "DesignOptions"))
    refstrat  <- which(colnames(x@StrataFrame)=="--")
    refstrat  <- if (length(refstrat)==0) 1L else refstrat[1]
    x@StrataFrame <- x@StrataFrame[refstrat]
    wtsum  <- sum(x@UnitWeights[complete.cases(x@StrataFrame)])
    if (is(x, "StratumWeightedDesignOptions"))
    {
        x@Sweights  <- x@Sweights[refstrat]
    } else {
        x  <- as(x, "StratumWeightedDesignOptions")
        ## The weights below are not meaningful -- even
        ## for the default harmonic weighting, there'd need
        ## to be a normalization factor -- but that's harmless
        ## because they'll have no effect on the shaping of
        ## covariates or stratum alignment provided by
        ## alignDesignsByStrata(). 
        x@Sweights <- setNames(
            list(data.frame(wtratio=rep(1, nlevels(x@StrataFrame[[1]])),
                                           row.names=levels(x@StrataFrame[[1]]))
                 ), names(x@StrataFrame)
        )
        }
    trans  <- if (is(center, "function")) center else NULL
    aligned  <-
        alignDesignsByStrata(colnames(x@StrataFrame)[1],
                             design=x,
                             post.align.transform=trans
                             )
    aligned_covs  <- aligned@Covariates
    if (scale)
        {
    scales  <- .colSums(aligned_covs^2,
                        nrow(aligned_covs), ncol(aligned_covs))
    scales  <- scales/wtsum
    scales[scales<.Machine$double.eps^.5]  <- 1
    sweep(aligned_covs, 2L, scales, "/", check.margin=FALSE)
        } else aligned_covs
    }

##' Align DesignOptions by Strata
##'
##' @param a_stratification name of a column of `design@strataFrame`/ element of `design@Sweights`
##' @param design DesignOptions
##' @param post.align.transform A post-align transform (cf \code{\link{balanceTest}})
##' @return CovsAlignedToADesign 
##' @keywords internal
##'
alignDesignsByStrata <- function(a_stratification, design, post.align.transform = NULL) {

    stopifnot(inherits(design, "StratumWeightedDesignOptions")) # defensive programming

    ss <- design@StrataFrame[, a_stratification]
    keep <- !is.na(ss)
    ss <- ss[keep]
    S <- SparseMMFromFactor(ss)
    Covs <- ifelse(design@NotMissing[keep, pmax(1L,design@NM.Covariates), drop = FALSE],
                   design@Covariates[keep, , drop = FALSE], 0)
    k.Covs <- ncol(Covs)
    NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
    vars  <- c(colnames(design@Covariates),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
    nmcols_Covs <- pmax(1L, design@NM.Covariates)
    origvars <- match(colnames(design@NotMissing), design@TermLabels, nomatch=0L)
    origvars <- c(design@OriginalVariables, origvars)

    UW  <- design@UnitWeights[keep]
    NM <- design@NotMissing[keep,,drop=FALSE]
    non_null_record_wts <- UW*NM[,1L,drop=TRUE]

    wtratio <- design@Sweights[[a_stratification]]$wtratio
    names(wtratio) <- rownames(design@Sweights[[a_stratification]])

    stopifnot(nlevels(ss)==1 ||
                  all(levels(ss)  %in% names(wtratio) ) )
    wtr.short <- wtratio[match(levels(ss),
                               names(design@Sweights[[a_stratification]]$wtratio),
                               nomatch=1L) # <-- this is to handle
                         ]                 # the unstratified case
    wtr <- wtr.short[ as.integer(ss) ]
    dim(wtr) <- NULL

    ## Figure weighted within-stratum stratum means, in order to 
    ## use them as imputation values. 
    Covs_stratmeans <- matrix(0, nrow(Covs), k.Covs)
    for (jj in 1L:k.Covs)
    {
        Covs_stratmeans[,jj] <- if (all(Covs[,jj]==Covs[1L,jj])) 0 else
            suppressWarnings( #warns if covar is linear in S or weights are all 0 w/in a stratum
                slm.wfit.csr( # see note in ./utils.R on why we use our own
                    S, Covs[,jj],               #`slm.wfit.csr` instead of `SparseM::slm.wfit`
                    weights=(UW*NM)[, nmcols_Covs[jj], drop = TRUE])$fitted
            )
    }
    ## Impute missing info to stratum mean
    Covs_w_touchups  <- Covs*NM[,nmcols_Covs, drop=FALSE] +
        Covs_stratmeans*(1-NM[,nmcols_Covs, drop=FALSE])

    if (!is.null(post.align.transform)) {
        ## Transform the columns of Covs_w_touchups using the post.align.trans
        Covs_aligned <- apply(Covs_w_touchups, 2,
                                 post.align.transform,
                                 non_null_record_wts #2nd arg to p.a.t.
                              )
        
    ## Ensure that post.align.trans wasn't something that changes the size of 
    ## Covs (e.g. mean). It would crash later anyway, but this is more informative
    if (is.null(dim(Covs_aligned)) || nrow(Covs_aligned) != nrow(Covs_w_touchups) ||
          ncol(Covs_aligned) != k.Covs) {
        stop("Invalid post.alignment.transform given")
      } else Covs_w_touchups  <- Covs_aligned
    }

    ## Align (recenter) cluster totals around their means within strata.
    ## Do this for the not-missing indicators as well as for the manifest variables.
    covars <- cbind(Covs_w_touchups, 0+NM[,NMcolperm])    
    covars <- suppressWarnings(
        slm_fit_csr(S, covars*non_null_record_wts)$residuals
    )
    colnames(covars) <- vars

      new("CovsAlignedToADesign",
          Covariates        = covars,
          Z=as.logical(design@Z[keep]),
          StrataMatrix=S,
          StrataWeightRatio = wtr, #as extracted from the design
          OriginalVariables = origvars,
          Cluster           = factor(design@Cluster[keep])
          )
}

# a helper function used by HB08 and HB08_2016 to see if problem is degenerate
# due to too many covariates
check_for_degenerate <- function(mat, n, s) {

   mat_rank <- attr(mat, "r")

  if (mat_rank >= (n - s)) {
    warning("Degrees of freedom exceeds units less number of strata, which can lead to a degenerate statistic. Try decreasing number of covariates tested.")
  }
}

##' @title Adjusted & combined differences as in Hansen & Bowers (2008)
##' @param alignedcovs A CovsAlignedToADesign object
##' @return list with components:
##' \describe{
##'   \item{z}{First item}
##'   \item{p}{Second item}
##'   \item{Msq}{Squared Mahalanobis distance of combined differences from origin}
##'   \item{DF}{degrees of freedom}
##'   \item{adj.diff.of.totals}{Vector of sum statistics z't - E(Z't), where t represents cluster totals of the product of the covariate with unit weights.  Hansen & Bowers (2008) refer to this as the adjusted difference vector, or d(z,x). }
##'   \item{tcov}{Matrix of null covariances of Z'x-tilde vector, as above.}
##' }
##' @references Hansen, B.B. and Bowers, J. (2008), ``Covariate
##'   Balance in Simple, Stratified and Clustered Comparative
##'   Studies,'' \emph{Statistical Science} \bold{23}.
##' @seealso \code{\link{balanceTest}}, \code{\link{alignDesignsByStrata}}
##' @importMethodsFrom SparseM diag
##' @keywords internal
HB08 <- function(alignedcovs) {
    zz <- as.numeric(alignedcovs@Z)
    S <- alignedcovs@StrataMatrix
    s_ <- ncol(S)

    Covs <- alignedcovs@Covariates
    n_  <- nrow(Covs)
    p_  <- ncol(Covs)

    ## by-stratum treatment and control counts
    n <- t(S) %*% S
    n.inv <- 1 / n
    n1 <- t(S) %*% zz
###    n0 <- t(S) %*% (1 - zz)
    n1_over_n <- S %*% n.inv %*% n1

    x_tilde <- Covs *#the sum statistic we're about to compute averages within-
        alignedcovs@StrataWeightRatio # stratum differences using stratum
    ## weights proportional to harmonic means of n_t's and n_c's.  If a different
    ## stratum weighting was indicated, it's shoehorned in here.  Whether
    ## or not that's so, weight normalization is also being factored in.
    ## In the default harmonic weighting, this weight ratio equals
    ## the reciprocal of the sum across strata of those harmonic weights.
    ## With another weighting, the weight ratios are all divided through
    ## by the sum of those other weights. 

    ssn <- sparseToVec(t(matrix(zz, ncol = 1) - n1_over_n) %*% x_tilde, column = FALSE)
    names(ssn) <- colnames(Covs)

    stratsizes <-  data.frame(n=diag(n), n1=sparseToVec(n1))
    stratsizes$n0  <- stratsizes[['n']] - stratsizes[['n1']]

    ## this next block first creates the n_ * p_^2 matrix
    ## of 2nd-order monomials in columns of x-tilde, then
    ## immediately sums each of these within each of s_ strata,
    ## resulting in a n_ * p_^2 matrix.
    xt_df  <- as.data.frame(x_tilde) # to get rep() to treat as a list
    xt_covar_stratwise  <- t(S) %*%
        ( as.matrix( as.data.frame(rep(xt_df, each=p_)) ) *
          as.matrix( as.data.frame(rep(xt_df, times=p_)) )
            )
    xt_covar_stratwise  <- as.matrix(xt_covar_stratwise) # s_ * (p_^2)
    xt_covar_stratwise  <- array(xt_covar_stratwise,
                                 dim=c(s_, p_, p_) # s_ * p_ * p_
                                 )
    xt_covar_stratwise  <-
        ifelse(stratsizes$n==1, 0, (stratsizes$n -1)^(-1) ) *
        xt_covar_stratwise # still s_ * p_ * p_
    ## now we have sample covariances, by stratum.

    ## Combine with factors equal to half the harmonic means of n1 and n0. 
    xt_c_s_scaled  <- xt_covar_stratwise *
        with(stratsizes, (1/n0 + 1/n1)^(-1) )
    tcov  <- apply(xt_c_s_scaled, 2:3, sum)

    ssvar <- diag(tcov)

    zero_variance  <- (ssvar <= .Machine$double.eps)
    zstat <- ifelse(zero_variance, NA_real_, ssn/sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)


    cov_minus_.5 <-
        XtX_pseudoinv_sqrt(mat=tcov[!zero_variance, !zero_variance, drop=FALSE],
                           mat.is.XtX = TRUE)

    check_for_degenerate(cov_minus_.5, n_, s_)

    mvz <- drop(crossprod(ssn[!zero_variance], cov_minus_.5))
    csq <- drop(crossprod(mvz))
    DF <- ncol(cov_minus_.5)

    list(z = zstat, p = p, Msq = csq , DF = DF,
         adj.diff.of.totals=ssn, tcov = tcov)
}

##' @title Hansen & Bowers (2008) inferentials 2016 [81e3ecf] version
##' @param alignedcovs A CovsAlignedToADesign object
##' @return list, as in \code{\link{HB08}}
##' @keywords internal
HB08_2016 <- function(alignedcovs) {
    zz <- as.numeric(alignedcovs@Z)
    S <- alignedcovs@StrataMatrix

    Covs <- alignedcovs@Covariates

    n <- t(S) %*% S
    n.inv <- 1 / n
    n1 <- t(S) %*% zz
    n0 <- t(S) %*% (1 - zz)

    # set up 1/(n-1)
    n_minus_1_inv <- n
    n_minus_1_inv@ra <- 1 / (n_minus_1_inv@ra - 1)
    n_minus_1_inv@ra <- ifelse(!is.finite(n_minus_1_inv@ra), # we can have strata of 1 unit,
                               0, n_minus_1_inv@ra) #  which causes a divide by zero error

    # product of {half the harmonic mean of n1, n0} with {1/(n-1)}
    dv <- sparseToVec(S %*% n_minus_1_inv %*% (n1 - n.inv %*% n1^2))

    x_tilde <- Covs *#the sum statistic we're about to compute averages within-
        alignedcovs@StrataWeightRatio # stratum differences using stratum
    ## weights proportional to harmonic means of n_t's and n_c's.  If a different
    ## stratum weighting was indicated, it's shoehorned in here. Whether
    ## or not that's so, weight normalization is also being factored in.
    ## In the default harmonic weighting, this weight ratio equals
    ## the reciprocal of the sum across strata of those harmonic weights.
    ## With another weighting, the weight ratios are all divided through
    ## by the sum of those other weights. 
    n1_over_n <- S %*% n.inv %*% n1
    ssn <- sparseToVec(t(matrix(zz, ncol = 1) - n1_over_n) %*% x_tilde, column = FALSE)
    names(ssn) <- colnames(Covs)

    scaled.x_tilde <- as.matrix(x_tilde * sqrt(dv))
    tcov <- crossprod(scaled.x_tilde)
    ssvar <- diag(tcov)

    zero_variance  <- (ssvar <= .Machine$double.eps)
    zstat <- ifelse(zero_variance, NA_real_, ssn/sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

    ## moving forward, we'll do without those sum statistics that have 0 null variation.
    x_tilde <- x_tilde[,zero_variance, drop=FALSE]
    scaled.x_tilde <- scaled.x_tilde[,ssvar > .Machine$double.eps, drop=FALSE]
    ssn  <- ssn[ssvar > .Machine$double.eps]

    cov_minus_.5 <- XtX_pseudoinv_sqrt(scaled.x_tilde)
    check_for_degenerate(cov_minus_.5, length(zz), ncol(S))
    
    mvz <- drop(crossprod(ssn, cov_minus_.5))
    csq <- drop(crossprod(mvz))
    DF <- ncol(cov_minus_.5)

    list(z = zstat, p = p, Msq = csq , DF = DF,
       adj.diff.of.totals=ssn, tcov = tcov)
}

##' Convert Matrix to vector
##'
##' @param s Matrix
##' @param column Column (TRUE) or row (FALSE)?
##' @return vector
##' @keywords internal
sparseToVec <- function(s, column = TRUE) {
  if (column) {
    as.matrix(s)[,1]
  } else {
    as.matrix(s)[1,]
  }
}
##' Identify vars recording not-missing (NM) info
##'
##' ID variables recording NM information, from
##' names and positions in the variable list.
##' Presumption is that the NM cols appear at the
##' end of the list of vars and are encased in
##' \sQuote{()}.  If something in the code changes
##' to make this assumption untrue, then this
##' helper is designed to err on the side of not
##' identifying other columns as NM cols.
##'
##' @param vnames character, variable names
##' @return character vector of names of NM vars, possibly of length 0
##' @author Hansen
##' @keywords internal
identify_NM_vars <- function(vnames) {
    k <- length(vnames)
    e <- nchar(vnames)
    inparens <- (substr(vnames,1L,1L)=="(" & substr(vnames,e,e)==")" )
    if (inparens[k] & any(!inparens) )
        {
            howmany <- which( cumsum( !inparens[k:1L] )==1 )[1] - 1
            vnames[(k+1-howmany):k]
            } else character(0)
    }
