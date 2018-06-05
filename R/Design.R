###############################################################################
# DesignMatrix Objects: covariates with term-specific missingness info & indexing
################################################################################

setClassUnion("Contrasts", c("list", "NULL"))
##' DesignMatrix S4 class
##'
##' If the Covariates matrix has an intercept, it will only be in the first column.
##'
##' More on NotMissing slot: It's matrix of numbers in [0,1].
##' First col is entirely \code{TRUE} or 1, like an intercept, unless corresponding
##' UnitWeight is 0, in which case it may also be 0 (see below). Subsequent cols
##' present only if there are missing covariate values, in which case these cols are
##' named for terms (of the original calling formula or data frame) that possess 
##' missing values.  Terms with the same missing data pattern are mapped to a single
##' column of this matrix.  If the DesignMatrix is representing elements, each column should
##' be all 1s and 0s, indicating which elements have non-missing values for the term
##' represented by that column.  If the DesignMatrix as a whole represents clusters,
##' then there can be fractional values, but that situation should only arise in the
##' DesignOptions class exension of this class, so it's documented there. 
##' 
##' @slot Covariates The numeric matrix that `model.matrix` would have returned.
##' @slot OriginalVariables look-up table associating Covariates columns with terms of the originating model formula
##' @slot term.labels labels of terms of the originating model formula
##' @slot contrasts Contrasts, a list of contrasts or NULL, as returned by `model.matrix.default`
##' @slot NotMissing Matrix of numbers in [0,1] with as many rows as the Covariates table but only one more col than there are distinct covariate missingness patterns (at least 1, nothing missing). First col is entirely T or 1, like an intercept.
##' @slot NM.Covariates integer look-up table mapping Covariates columns to columns of NotMissing.  (If nothing missing for that column, this is 0.)
##' @slot NM.terms integer look-up table mapping term labels to columns of NotMissing (0 means nothing missing in that column)
##' @keywords internal
##' 
setClass("DesignMatrix",
         slots=c(Covariates="matrix",
                  OriginalVariables="integer",
                  TermLabels="character",
                  Contrasts="Contrasts",
                  NotMissing="matrix",
                  NM.Covariates="integer",
                  NM.terms="integer" )
         )

#' @export
as.matrix.DesignMatrix <- function(x, ...)
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
##' @param data Data frame (as in `model.matrix`)
##' @param remove.intercept logical
##' @param ... passed to `model.matrix.default` (and further)
##' @return DesignMatrix instance of an S4 class that enriches model matrices with missing data info
##' @author Ben B Hansen
##' @keywords internal
##' 
design_matrix <- function(object, data = environment(object), remove.intercept=TRUE, ...) {
  # mf <- model.frame(object, data, na.action = na.pass)
  tms <- terms(object)
  term.labels <- attr(tms, "term.labels")

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
    
  ccs.by.term <- c(list('_non-null record_'=!null.record),
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
      ## _non-null record_ and for each provided term. Entries are
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
    
  new("DesignMatrix",
      Covariates=covariates,
      OriginalVariables=assign,
      TermLabels=term.labels,
      Contrasts=contrasts,
      NotMissing=notmissing,
      NM.Covariates=nm.covs,
      NM.terms=nm.terms )
}


################################################################################
# DesignOptions Objects: communicate cluster, strata, treatment, and covariate information
################################################################################

#' DesignOptions S4 class
#'
#' Extends the DesignMatrix class
#'
##' If the DesignOptions represents clusters of elements, as when it was created
##' by aggregating another DesignOptions or DesignMatrix object, then its Covariates
##' and NotMissing slots are populated with (weighted) averages, not totals.  E.g.,
##' NotMissing columns consist of weighted averages of element-wise non-missingness indicators
##' over clusters, with weights given by (the element-level precursor to) the UnitWeights
##' vector.  As otherwise, columns of the NotMissing matrix represent terms
##' from a model formula, rather than columns the terms may have expanded to.  

#' @slot Z Logical indicating treatment assignment
#' @slot StrataMatrices This is a list of sparse matrices, each with n rows and s columns, with 1 if the unit is in that stratification
#' @slot StrataFrame Factors indicating strata (not the sparse matrices, as we use them in the weighting function)
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @slot UnitWeights vector of weights associated w/ rows of the DesignMatrix
#' @keywords internal
#' 
setClass("DesignOptions",
         representation = list(
           Z                 = "logical",
           StrataMatrices    = "list", 
           StrataFrame       = "data.frame",  
             Cluster           = "factor",
             UnitWeights = "numeric"),
         contains = "DesignMatrix"
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

    uweights <- as.vector(model.weights(data))
    if (is.null(uweights))
        stop("makeDesigns() expects its data arg to be a model frame containing weights")
    stopifnot(is.numeric(uweights), all(!is.na(uweights)), all(uweights>=0))

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
  # Following resolution to #86 in [master ad6ed6a], we have to indicate specifically
  # that `cluster` and `strata` are to be found in the survival package.
  str.vnames.safe <- gsub('(?<!:)cluster\\(', 'survival::cluster\\(', str.vnames, perl=TRUE)
  str.vnames.safe <- gsub('(?<!:)strata\\(', 'survival::strata\\(', str.vnames.safe, perl=TRUE)     
  # The purposes of the regexp lookbehinds (`(?<!:)`) above are to avoid overwriting
  # "survival::cluster(" with "survival::survival::cluster(", and also to avoid
  # overruling users who prefer to get their `cluster()` or `strata()` from elsewhere 
  # than the survival package. 

  str.fmla <- formula(paste0("factor(", treatment.name, ")", " ~ ", paste0(collapse = "+", c(1, str.vnames.safe))),
                      env=environment(fmla))
  str.tms  <- terms(str.fmla, data = data,
                    specials = c("survival::cluster", "survival::strata"))
  str.data <- model.frame(str.tms, data = data, na.action = na.pass, drop.unused.levels=TRUE)


  ## check that strata and clusters have the proper relationships with treatment assignment
  treatmentCol <- colnames(str.data)[attr(str.tms, "response")]
  clusterCol <- colnames(str.data)[attr(str.tms, "specials")$`survival::cluster`]
  strataCols <- colnames(str.data)[attr(str.tms, "specials")$`survival::strata`]

  if (includeUnstratified) {
    str.data$Unstrat <- 1
    strataCols <- c(strataCols, "Unstrat")
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
    str.data[ str.data[, s]%in%colnames(tbl)[!isGood], s] <- NA

    if (!(all(isGood))) {
        warning("Dropped ", sum(!isGood), " levels of ", s, " which did not include both treated and control units")
###      warning("In ", s, ", dropping the following stratum levels (which do not include both treated and control units):\n",
###           paste(collapse = ", ", colnames(tbl)[!isGood]))
  }
  }

  ## OK! data looks good. Let's proceed to make the design object with covariate data

  data.fmla <- update(ts, paste("~", paste0(collapse = " - ", c(".", str.vnames))))
  data.data <- model.frame(data.fmla, data, na.action = na.pass) #

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

    desm <- design_matrix(terms(data.data), data.data, remove.intercept=TRUE, contrasts.arg = clist)

  if (length(clusterCol) > 0) {
    Cluster <- str.data[, clusterCol]
  } else {
    Cluster <- 1:(dim(str.data)[1])
  }

  Z <- str.data[, treatmentCol]
  tmp <- str.data[, strataCols, drop = FALSE]
  colnames(tmp) <- gsub(colnames(tmp), pattern = "survival::strata\\((.*)\\)", replacement = "\\1")
  strata.frame <- data.frame(lapply(tmp, factor), check.names = FALSE)
  strata.mats  <- lapply(strata.frame, function(s) { SparseMMFromFactor(s) })
  names(strata.mats) <- colnames(strata.frame)


  return(new("DesignOptions",
             Z                 = as.logical(as.numeric(Z) - 1), #b/c it was built as factor
             StrataMatrices    = strata.mats,
             StrataFrame       = strata.frame,
             Cluster           = factor(Cluster),
             UnitWeights = uweights,
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
##' culled or inferred from originating \link{\code{balanceTest}} call or as aggregated up
##' from those unit weights).  Returns a 
##' weighting factor to be associated with each stratum, this factor determining the stratum 
##' weight by being multiplied by mean of unit weights over clusters in that stratum.
##'
##' Specifically, the function's value is a data frame of two variables,
##' \code{sweights}  and \code{wtratio}, with rows representing strata.
##' The \code{sweights} vector represents internally
##' calculated or user-provided \code{stratum.weights}, scaled so that
##' its sum over each of the strata is 1.  Other than this scaling
##' it's w_b of Hansen & Bowers (2008). \code{wtratio} is the ratio of
##' \code{sweights} to the product of the harmonic
##' mean of n_{tb} and n_{cb}, the number of treatment and control
##' clusters in stratum b, with the mean of the weights associated with
##' each of these clusters,  i.e. the m-bar_b of Hansen & Bowers
##' (2008). This comparison of \code{sweights} to the product of h_b
##' and m-bar_b is expected downstream in
##' \code{AlignedToInferentials} (in its internal calculations
##' involving \sQuote{\code{wtr}}).
##' 
##' (Developer note: One might simplify by only returning the sweights,
##' not also the wtratio's.  [Or perhaps conversely.]  sweights are used in descriptives
##' calculations.  wtratio's are used downstream in 
##' alignedToInferentials, but the material in the harmonic means calcs is already 
##' being assembled there for other reasons.  An improvement for another day...)
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

    if (nlevels(stratifier) == 1) {
      return(data.frame(sweights=1, wtratio=1, row.names='1'))
    }

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
      wtratio <- rep(1, length(swts))
    } else {
        hwts <- harmonic_times_mean_weight(
            data.frame(Tx.grp = design@Z,
                       stratum.code=stratifier,
                       unit.weights=design@UnitWeights,
                       check.names = FALSE))
      wtratio <- swts/hwts
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
##' @param design A DesignOptions object
##' @param covariate.scaling Scale estimates for covs, to use instead of internally calculated pooled SDs
##' @return Descriptives
##' @keywords internal
##' 
designToDescriptives <- function(design, covariate.scaling = NULL) {
  stopifnot(inherits(design, "DesignOptions")) # defensive programming
  if (!is.null(covariate.scaling)) warning("Non-null 'covariate.scaling' currently being ignored")
  covars <- ifelse(is.na(design@Covariates), 0, design@Covariates)

  ## Tack NM cols onto covars, but with intercept col listed last
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  covars <- cbind(covars, design@NotMissing[, NMcolperm])
  vars <- c(colnames(design@Covariates),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
  colnames(covars)   <-  vars
  covars.nmcols <- c(pmax(1L, design@NM.Covariates), rep(1L, k.NM ) )
  stratifications <- colnames(design@StrataFrame)

  Uweights <- design@UnitWeights * design@NotMissing

  ans <- array(NA,
               dim = c(length(vars), 5, length(stratifications)),
               dimnames = list(
                   "vars" = vars,
                   "stat" = c("Control", "Treatment", "std.diff", "adj.diff", "pooled.sd"),
                   "strata" = stratifications))
  for (s in stratifications) {

    S <- design@StrataMatrices[[s]]
    if (ncol(S)!=nlevels(design@StrataFrame[[s]]))
        stop(paste("Levels of StrataFrame don't match StratMatrices colnames, stratification", s))
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

    S.missing.0 <- as.matrix((t(ZZ) %*% Uweights)) == 0
    S.missing.1 <- as.matrix((t(WW) %*% Uweights)) == 0
    S.has.both  <- !(S.missing.0 | S.missing.1)
    use.units   <- S %*% S.has.both * Uweights
    use.units <- as.matrix(use.units)

    X.use  <- covars * use.units[, covars.nmcols]
    X2.use <- covars^2 * use.units[, covars.nmcols]

    n1 <- t(use.units) %*% Z
    n1 <- n1[covars.nmcols, ]
    n0 <- t(use.units) %*% (1 - Z)
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
    wtsum.tx <-  t(use.units * tx.wts) %*% Z
    wtsum.tx <- wtsum.tx[covars.nmcols, ]
    treated.avg <- t(X.use *tx.wts) %*% Z / wtsum.tx

    ## ratio estimates of means over similarly defined domains within the control group
    wtsum.ctl <- t(use.units * ctl.wts) %*% (1 - Z)
    wtsum.ctl <- wtsum.ctl[covars.nmcols, ]
    control.avg <- t(X.use * ctl.wts) %*% (1 - Z) / wtsum.ctl

    var.1 <- (t(X2.use *tx.wts) %*% Z - wtsum.tx * treated.avg^2) / wtsum.tx
    var.1 <- var.1 * ifelse(n1>1, n1/(n1 - 1), 0)
    var.0 <- (t(X2.use * ctl.wts) %*% (1 - Z) - wtsum.ctl * control.avg^2) / wtsum.ctl
    var.0 <- var.0* ifelse(n0>1, n0/(n0 - 1), 0)

    pooled <- sqrt((var.1 + var.0) / 2)

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
##' @param design DesignOptions
##' @return another DesignOptions representing the clusters
##' @keywords internal
##' 
aggregateDesigns <- function(design) {
  n.clusters <- nlevels(design@Cluster)

  if (n.clusters == length(design@Cluster)) {
    return(design)
  }

  dupes <- duplicated(design@Cluster)
  Z <- design@Z[!dupes]
  StrataFrame <- design@StrataFrame[!dupes,, drop = FALSE]
  Cluster <- design@Cluster[!dupes]

  C <- SparseMMFromFactor(design@Cluster)

  unit.weights <- as.matrix(t(C) %*% as.matrix(design@UnitWeights))
  dim(unit.weights) <- NULL
  
  Uweights.tall <- design@UnitWeights * design@NotMissing
  Covariates <- as.matrix(t(C) %*% ifelse(Uweights.tall[,pmax(1L,design@NM.Covariates), drop=FALSE],
                                          design@Covariates *
                                              Uweights.tall[,pmax(1L,design@NM.Covariates), drop=FALSE],
                                          0)
                          )
  Uweights <- as.matrix(t(C) %*% Uweights.tall)
  Covariates <- ifelse(Uweights[,pmax(1L,design@NM.Covariates), drop=FALSE] > 0,
                       Covariates/Uweights[,pmax(1L,design@NM.Covariates), drop=FALSE],
                       0)
  NotMissing <- ifelse(matrix(unit.weights>0, nrow(Uweights), ncol(Uweights)),
                       Uweights/unit.weights, 0)
  colnames(NotMissing) <- colnames(design@NotMissing)
  
  StrataMatrices <- lapply(design@StrataMatrices, function(S) {
    tmp <- t(C) %*% S
    tmp@ra <- rep(1, length(tmp@ra))
    return(tmp)
  })

  # colnames(Covariates) <- c("cluster.size", colnames(design@Covariates))
  colnames(Covariates)   <- colnames(design@Covariates)

  new("DesignOptions",
      Z = Z,
      StrataMatrices = StrataMatrices,
      StrataFrame = StrataFrame,
      Cluster = Cluster,
      UnitWeights = unit.weights,
      NotMissing = NotMissing,
      Covariates = as.matrix(Covariates),
      OriginalVariables = design@OriginalVariables,
      TermLabels=design@TermLabels,
      Contrasts=design@Contrasts,
      NM.Covariates=design@NM.Covariates,
      NM.terms=design@NM.terms)
}

#' AlignedCovs S4 class
#'
#' A class for representing covariate matrices after alignment within stratum,
#' for a (single) given stratifying factor.  There can also be a clustering variable,
#' assumed to be nested within the stratifying variable. Extends DesignMatrix class.
#'
#' Unlike DesignOptions this class has no UnitWeights slot.  Rather, the NotMissing
#' slot now records products of unit weights and non-missingness (more specifically
#' cluster means of non-missingness averaged with element weights, if provided). 
#' The StrataWeightRatio slot has an entry for each unit, representing ratio of
#' specified stratum weight to the product of h_b (the harmonic mean of n_{tb} and
#' n_{cb}, the counts of treatment and control clusters in stratum b) with bar-w_b,
#' (the arithmetic mean of aggregated cluster weights within that stratum).
#' 
#' @slot Z Logical indicating treatment assignment
#' @slot StrataMatrix A sparse matrix with n rows and s columns, with 1 if the unit is in that stratification
#' @slot StrataFactor Factor indicating strata
#' @slot StrataWeightRatio For each unit, ratio of stratum weight to h_b; but see Details.
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @slot OriginalVariables Look up table associating Covariates cols to terms in the calling formula, as in DesignMatrix
#' @slot Covariates Numeric matrix, as in DesignMatrix, except here we presume columns to have been aligned (stratum-centered)
#' @slot NotMissing matrix of unit weights, normalized within stratum to have mean 1. (As otherwise, NAs are 0s)
setClass("AlignedCovs",
         slots = list(
             Z                 = "logical",
             StrataMatrix    = "matrix.csr", 
             StrataFactor       = "factor",
             StrataWeightRatio = "numeric",
             Cluster = "factor"
         ),
         contains="DesignMatrix")



##' Align DesignOptions by Strata
##'
##' @param design DesignOptions
##' @param post.align.transform A post-align transform
##' @return list List of `AlignedCovs` objects
##' @keywords internal
##' 
alignDesignsByStrata <- function(design, post.align.transform = NULL) {

  stopifnot(inherits(design, "StratumWeightedDesignOptions")) # defensive programming

  vars   <- colnames(design@Covariates)
  strata <- names(design@StrataMatrices)

  Ewts  <- design@UnitWeights * design@NotMissing
  Covs <- ifelse(design@NotMissing[, pmax(1L,design@NM.Covariates), drop = FALSE],
                 design@Covariates[, , drop = FALSE], 0)
  k.Covs <- ncol(Covs)
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  Covs <- cbind(Covs, 0+design@NotMissing[,NMcolperm])
  vars  <- c(colnames(design@Covariates),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
  covars.nmcols <- c(design@NM.Covariates, rep(1L, k.NM ) )
  origvars <- match(colnames(design@NotMissing), design@TermLabels, nomatch=0L)
  origvars <- c(design@OriginalVariables, origvars)

  # we can't return an array because different stratifications will have varying numbers
  # of strata levels. A list is more flexible here, but less structured.
  f <- function(s){
    ss <- design@StrataFrame[, s]
    keep <- !is.na(ss)
    ss <- ss[keep]
    S <- SparseMMFromFactor(ss)

    ewts <- Ewts[keep,,drop=FALSE]
    NM <- design@NotMissing[keep,,drop=FALSE]
    NMCovs <- design@NM.Covariates
    covars <- Covs[keep,,drop=FALSE]

    wtratio <- design@Sweights[[s]]$wtratio
    names(wtratio) <- rownames(design@Sweights[[s]])

    stopifnot(nlevels(ss)==1 ||
                  all(levels(ss)  %in% names(wtratio) ) )
    wtr.short <- wtratio[match(levels(ss),
                               names(design@Sweights[[s]]$wtratio),
                               nomatch=1L) # <-- this is to handle 
                         ]                 # the unstratified case
    wtr <- wtr.short[ as.integer(ss) ]
    dim(wtr) <- NULL

    # align weighted observations within stratum by subtracting weighted stratum means
    covars.Sctr <- covars
    for (jj in 1L:ncol(covars))
    {
        covars.Sctr[,jj] <- if (all(covars[,jj]==covars[1L,jj])) 0 else 
            suppressWarnings( #throws singularity warning if covar is linear in S
                slm.wfit.csr( # see note in ./utils.R on why we use our own 
                    S, covars[,jj],               #`slm.wfit.csr` instead of `SparseM::slm.wfit`
                    weights=ewts[, max(1L, covars.nmcols[jj]), drop = TRUE])$residuals
            )
        ## A value that was missing might have received an odd residual.
        ## Although such values don't themselves contribute anything, they'll affect a
        ## post alignment transformation such as `rank`.  So, per #47 we set them to 0 (the stratum mean).
        if (jj<= k.Covs)
            covars.Sctr[ !NM[, NMCovs[jj] ], # picks out rows w/ missing observations
                   jj] <- 0
    }

    if (!is.null(post.align.transform)) {
        ## Transform the columns of covars.Sctr using the function in post.align.trans
        ## only do so for actual covariates, however, not missingness weights
      covars.Sctr.new <- apply(covars.Sctr[,1:k.Covs, drop=FALSE], 2, post.align.transform)

      # Ensure that post.align.trans wasn't something that changes the size of covars.Sctr (e.g. mean).
      # It would crash later anyway, but this is more informative
      if (is.null(dim(covars.Sctr.new)) || nrow(covars.Sctr) != nrow(covars.Sctr.new) ||
          ncol(covars.Sctr.new) != k.Covs) {
        stop("Invalid post.alignment.transform given")
      }
      ## The post alignment transform may have disrupted the stratum alignment.  So, recenter on stratum means
        covars.Sctr[,1L:k.Covs] <- suppressWarnings(
            slm.wfit.csr(S, covars.Sctr.new[,1L:k.Covs, drop=FALSE],
                         weights=ewts[, 1L, drop = TRUE])$residuals
        )
    }
    colnames(covars.Sctr) <- vars
      
      new("AlignedCovs",
          Z=as.logical(design@Z[keep]),
          StrataMatrix=S,
          StrataFactor=ss,
          StrataWeightRatio = wtr, #as extracted from the design
          Cluster           = factor(design@Cluster[keep]),
          OriginalVariables = origvars,
          Covariates        = covars.Sctr,
          NotMissing          = ewts,
          TermLabels=design@TermLabels,
          Contrasts=design@Contrasts,
          NM.Covariates=covars.nmcols,
          NM.terms=design@NM.terms)
}
  sapply(strata, f, simplify = FALSE, USE.NAMES = TRUE)
}

# I'd prefer to have a better API here, but right now, just trying to get compatability with old xBalance set up
# e.g. something that is a list of strata with a given structure, rather than just a list.
##' Align to Inferentials
##'
##' @param alignedcovs An AlignedCovs object
##' @return list
alignedToInferentials <- function(alignedcovs) {
    zz <- as.numeric(alignedcovs@Z)
    S <- alignedcovs@StrataMatrix
    wtr <- alignedcovs@StrataWeightRatio

    ## we need to map the covariates to their columns in the NotMissing matrix.
    ## if a column has no missing at all, we indicate that with a zero
    ## but this would otherwise cause the column to get dropped.
    ## Instead, we map it to the first column of NM.Covariates, which we know is all 1s.
    mapping <- alignedcovs@NM.Covariates
    mapping[mapping == 0] <- 1
    Uweights <- alignedcovs@NotMissing[,mapping]
    Covs <- alignedcovs@Covariates
    
    n <- t(S) %*% S
    n.inv <- 1 / n
    n1 <- t(S) %*% zz
    n0 <- t(S) %*% (1 - zz)

    # dv is sample variance of treatment by stratum
    # set up 1/(n-1)
    tmp <- n
    tmp@ra <- 1 / (tmp@ra - 1)
    dv <- sparseToVec(S %*% tmp %*% (n1 - n.inv %*% n1^2)) 
    
    tmat <- Covs * Uweights * #the sum statistic we're about to compute corresponds 
        wtr # to averaging within-stratum difference w/ stratum weights proportional to
    ## harmonic means of n_ts and n_cs.  To override w/user-designated weights, we
    ## factor in wtr, the "weight ratio" as previously reconstructed.

    ZtH <- S %*% n.inv %*% n1
    ssn <- sparseToVec(t(matrix(zz, ncol = 1) - ZtH) %*% tmat, column = FALSE)
    scaled.tmat <- as.matrix(tmat * sqrt(dv))
    tcov <- crossprod(scaled.tmat)
    ssvar <- diag(tcov)    

    ## The next few calcs put components of 
    ## z-statistic onto an interpretable scale      
    ##  wtsum is the sum across strata of twice the harmonic mean of n1, n0 - we should rename it
    wtsum <- sum((n.inv %*% (n1 * n0))@ra) # (the ra slot is where SparseM keeps the non-zero values)
    post.diffs <- ssn / wtsum
    tcov <- tcov *(1 / wtsum^2)
    
    zstat <- ifelse(ssvar <= .Machine$double.eps, NA_real_, ssn/sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

    ## moving forward, we'll do without those sum statistics that have 0 null variation.
    tmat <- tmat[,ssvar > .Machine$double.eps, drop=FALSE]
    scaled.tmat <- scaled.tmat[,ssvar > .Machine$double.eps, drop=FALSE]
    pst.svd <- try(svd(scaled.tmat, nu=0))

  if (inherits(pst.svd, 'try-error')) {
    pst.svd <- propack.svd(scaled.tmat)
  }

  Positive <- pst.svd$d > max(sqrt(.Machine$double.eps) * pst.svd$d[1], 0)
  Positive[is.na(Positive)] <- FALSE 

  if (all(Positive)) { ## is this faster? { ytl <- sweep(pst.svd$v,2,1/pst.svd$d,"*") }
    ytl <- pst.svd$v *
      matrix(1/pst.svd$d, nrow = dim(tmat)[2], ncol = length(pst.svd$d), byrow = T)
  } else if (!any(Positive)) {
    ytl <- array(0, dim(tmat)[2:1] )
  } else  {
    ytl <- pst.svd$v[, Positive, drop = FALSE] *
      matrix(1/pst.svd$d[Positive], ncol = sum(Positive), nrow = dim(tmat)[2], byrow = TRUE)
  }

  mvz <- drop(crossprod(zz, tmat) %*% ytl)

  csq <- drop(crossprod(mvz))
  DF <- sum(Positive)

  list(z = zstat, p = p, csq = csq , DF = DF, 
       adj.mean.diffs=post.diffs, tcov = tcov)
}

##' Convert Matrix to vector
##'
##' @param s Matrix
##' @param column Column (TRUE) or row (FALSE)?
##' @return vector
sparseToVec <- function(s, column = TRUE) {
  if (column) {
    as.matrix(s)[,1]
  } else {
    as.matrix(s)[1,]
  }
}

