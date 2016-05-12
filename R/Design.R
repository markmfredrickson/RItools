################################################################################
# DesignMatrix Objects: covariates with term-specific missingness info & indexing
################################################################################

setClassUnion("Contrasts", c("list", "NULL"))
##' DesignMatrix S4 class
##'
##' If the Covariates matrix has an intercept, it will only be in the first column.
##' @slot Covariates The numeric matrix that `model.matrix` would have returned.
##' @slot OriginalVariables look-up table associating Covariates columns with terms of the originating model formula
##' @slot term.labels labels of terms of the originating model formula
##' @slot contrasts Contrasts, a list of contrasts or NULL, as returned by `model.matrix.default`
##' @slot NotMissing Matrix of numbers in [0,1] with as many rows as Covariates but only one more col than there are distinct covariate missingness patterns (at least 1, nothing missing). First col is entirely T or 1, like an intercept.
##' @slot NM.Covariates integer look-up table mapping Covariates columns to columns of NotMissing.  (If nothing missing for that column, this is 1.)
##' @slot NM.terms integer look-up table mapping term labels to  columns of NotMissing
setClass("DesignMatrix",
         slots=c(Covariates="matrix",
                  OriginalVariables="integer",
                  TermLabels="character",
                  Contrasts="Contrasts",
                  NotMissing="matrix",
                  NM.Covariates="integer",
                  NM.terms="integer" )
         )

##' @export
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
design_matrix <- function(object, data=environment(object), remove.intercept=TRUE, ...)
    {
        mf <- model.frame(object, data, na.action=na.pass)
        tms <- attr(mf, "terms")
        term.labels <- attr(tms, "term.labels")

        covariates <- model.matrix(object=object, data=mf, ...) 
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
        

        cols.by.term <- lapply(1L:length(term.labels),
                               function(whichterm) (1:ncol(covariates))[assign==whichterm])
        ccs.by.term <-
            lapply(cols.by.term, function(whichcols) {
                complete.cases(covariates[,whichcols])
            })
        terms.with.missings <- !sapply(ccs.by.term, all)

        nm.covs <- integer(ncol(covariates))
        nm.terms <- integer(length(term.labels))

        if (any(terms.with.missings)) {
            nmdf <- as.data.frame(ccs.by.term[terms.with.missings])
            colnames(nmdf) <-  term.labels[terms.with.missings]
            nmdf.char <- sapply(nmdf, function(x) paste(as.integer(x), collapse="."))
            nmdf.dupes <- duplicated.default(nmdf.char, fromLast=FALSE)
            if (any(nmdf.dupes))
                {
                    nmdf <- nmdf[!nmdf.dupes]
                    nm.terms[terms.with.missings][!nmdf.dupes] <- 1L:ncol(nmdf)
                    nm.terms[terms.with.missings][nmdf.dupes] <-
                        match(nmdf.char[nmdf.dupes], nmdf.char[!nmdf.dupes])
                } else nm.terms[terms.with.missings] <- 1L:ncol(nmdf)
            nm.covs[assign>0] <- nm.terms[assign[assign>0]]
            notmissing <- as.matrix(nmdf)
        } else {
            notmissing <- matrix(FALSE, nrow(covariates), 0)
        }

        ## add in a 1st column of 1s, to ease bookkeeping later on
        notmissing <- cbind(matrix(TRUE, nrow(covariates), 1), notmissing)
        colnames(notmissing)[1] <- "element weight"
        nm.covs <- nm.covs + 1L
        nm.terms <- nm.terms + 1L

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
#' @slot Z Logical indicating treatment assignment
#' @slot StrataMatrices This is a list of sparse matrices, each with n rows and s columns, with 1 if the unit is in that stratification
#' @slot StrataFrame Factors indicating strata (not the sparse matrices, as we use them in the weighting function)
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @slot NotMissing Matrix of numbers in [0,1], cluster means of element weights with missing values
#' contributing 0s.  First col is entirely T or 1, like an intercept, unless corresponding
#' ElementWeight is 0, in which case it may also be 0. Subsequent cols present only if
#' there are missing covariate values, in which case these cols are named for terms that
#' possess missing values.  Terms with the same missing data pattern are mapped to a single
#' 
setClass("DesignOptions",
         representation = list(
           Z                 = "logical",
           StrataMatrices    = "list", 
           StrataFrame       = "data.frame",  
             Cluster           = "factor",
             ElementWeights = "numeric"),
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
makeDesigns <- function(fmla, data) {

    eweights <- as.vector(model.weights(data))
    if (is.null(eweights))
        stop("makeDesigns() expects its data arg to be a model frame containing weights")
    stopifnot(is.numeric(eweights), all(!is.na(eweights)), all(eweights>=0))

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
  str.idx <- c(attr(ts, "specials")$cluster, attr(ts, "specials")$strata)

  str.fmla <- formula(paste0("factor(", treatment.name, ")", " ~ ", paste0(collapse = "+", c(1, vnames[str.idx]))))
  str.tms  <- terms(str.fmla, data = data, specials = c("cluster", "strata"))
  str.data <- model.frame(str.tms, data = data, na.action = na.pass, drop.unused.levels=TRUE)


  ## check that strata and clusters have the proper relationships with treatment assignment
  treatmentCol <- colnames(str.data)[attr(str.tms, "response")]
  clusterCol <- colnames(str.data)[attr(str.tms, "specials")$cluster]
  strataCols <- colnames(str.data)[attr(str.tms, "specials")$strata]

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

  data.fmla <- update(ts, paste("~", paste0(collapse = " - ", c(".", vnames[str.idx]))))
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
  colnames(tmp) <- gsub(colnames(tmp), pattern = "strata\\((.*)\\)", replacement = "\\1")
  strata.frame <- data.frame(lapply(tmp, factor), check.names = FALSE)
  strata.mats  <- lapply(strata.frame, function(s) { SparseMMFromFactor(s) })
  names(strata.mats) <- colnames(strata.frame)


  return(new("DesignOptions",
             Z                 = as.logical(as.numeric(Z) - 1), #b/c it was built as factor
             StrataMatrices    = strata.mats,
             StrataFrame       = strata.frame,
             Cluster           = factor(Cluster),
             ElementWeights = eweights,
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
##' (Developer note: There's no real reason not to simplify by only returning the sweights,
##' not also the wtratio's.  Weight ratios are used downstream in alignedToInferentials, but
##' the material in the harmonic means calcs is already being assembled there for other
##' reasons.  An improvement for another day...)
##' @param design DesignOptions
##' @param stratum.weights Stratum weights
##' @return list with two entries, \code{sweights}  and \code{wtratio}.
##' each is in turn a list with an entry for each stratum

DesignWeights <- function(design, stratum.weights = harmonic) {
  stopifnot(inherits(design, "DesignOptions"))

  n.strata <- dim(design@StrataFrame)[2]
  strata.names <- colnames(design@StrataFrame)

  if (is.function(stratum.weights)) {
    swt.ls <- rep(list(stratum.weights), n.strata)
    names(swt.ls) <- strata.names
  }


  if (!is.list(stratum.weights) & !is.function(stratum.weights) & n.strata > 1)
    stop("stratum weights must be specified for each stratifying factor")

  if (!is.list(stratum.weights) & !is.function(stratum.weights)) {
    swt.ls <- list(stratum.weights)
    names(swt.ls) <- strata.names
  }

  if (is.list(stratum.weights) & !all(strata.names %in% c("Unstrat", names(stratum.weights))))
    stop("list stratum.weights must have entry names matching those of stratifying factors")

  if (is.list(stratum.weights) & !is.function(stratum.weights)) {

      if (("Unstrat" %in% strata.names) && !("Unstrat" %in% names(stratum.weights)) )
          stratum.weights <- c(stratum.weights, list("Unstrat"=NULL))

      swt.ls <- stratum.weights[strata.names]
  }


  ### change names here!

  wtlist <- list()
  for (nn in names(swt.ls)) {

      if (nlevels(factor(design@StrataFrame[[nn]]))==1)
          {
              wtlist[[nn]] <- data.frame(sweights=1, wtratio=1, row.names='1')
              next
          }

    if (is.function(swt.ls[[nn]])) {
      sweights <-
        do.call(swt.ls[[nn]],
                args = list(data =
                    data.frame(Tx.grp = design@Z,
                               stratum.code = factor(design@StrataFrame[[nn]]),
                               design@Covariates,
                               check.names = FALSE)),
                envir=parent.frame())
      Eweight.stratum.means <- tapply(design@ElementWeights, design@StrataFrame[[nn]], mean)
      sweights <- Eweight.stratum.means * sweights
    } else {
      if (!is.numeric(swt.ls[[nn]]))
        stop("stratum.weights must be an expression or numeric vector")

      if (is.null(names(swt.ls[[nn]])))
        stop ("if stratum.weights is a vector, must have names")

      if (!(all(levels(factor(design@StrataFrame[[nn]])) %in% names(swt.ls[[nn]])) ))
        stop("if stratum.weights is a vector, must have a name for each stratum")

      sweights <- swt.ls[[nn]][levels(factor(design@StrataFrame[[nn]]))]
    }

    if (all(is.na(sweights)))
      stop(paste("All stratum weights NA (strat.",nn,")."))
    if (any(is.na(sweights))) {
      sweights[is.na(sweights)] <- 0
      warning(paste("NAs in stratum weights (",nn," strat.); to be interpreted as 0s.", sep=""))
    }
    if (any(sweights<0))
      stop("stratum weights must be nonnegative")

      sweights <- sweights/sum(sweights, na.rm=TRUE)

    if (identical(harmonic, swt.ls[[nn]])) {
      hwts <- sweights
    } else {
      hwts <- harmonic(data.frame(Tx.grp = design@Z,
                                  stratum.code=factor(design@StrataFrame[[nn]]),
                                  check.names = FALSE))
    }
    hwts <- hwts/sum(hwts, na.rm=TRUE)

      wtratio <- sweights/hwts
#    wtratio <- unsplit(sweights/hwts, design@StrataFrame[[nn]], drop=TRUE)
#    wtratio[is.na(wtratio)] <- 0
    wtlist[[nn]] <- data.frame(sweights=sweights,wtratio=wtratio)
    NULL
  }
  wtlist
}

##' Generate Descriptives
##'
##' Use a design object to generate descriptive statistics that ignore clustering. 
##' Stratum weights are respected.
##' @param design A DesignOptions object
##' @param covariate.scaling Scale estimates for covs, to use instead of internally calculated pooled SDs
##' @return Descriptives
designToDescriptives <- function(design, covariate.scaling = NULL) {
  stopifnot(inherits(design, "DesignOptions")) # defensive programming
  if (!is.null(covariate.scaling)) warning("Non-null 'covariate.scaling' currently being ignored")
  covars <- ifelse(is.na(design@Covariates), 0, design@Covariates)

  ## Tack NM cols onto covars, but with intercept col listed last
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  covars <- cbind(covars, design@NotMissing[, NMcolperm])
  vars <- c(colnames(design@Covariates),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
  colnames(covars)   <-  vars
  covars.nmcols <- c(design@NM.Covariates, rep(1L, k.NM ) )
  stratifications <- colnames(design@StrataFrame)

  Eweights <- design@ElementWeights * design@NotMissing
#  Eweights <- Eweights


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
            Swts <- rep(1, length(stratlevs))
            names(Swts) <- stratlevs
                }
    
    Z <- as.numeric(design@Z)
    ZZ <- S * Z
    WW <- S * (1 - Z)

    S.missing.0 <- as.matrix((t(ZZ) %*% Eweights)) == 0
    S.missing.1 <- as.matrix((t(WW) %*% Eweights)) == 0
    S.has.both  <- !(S.missing.0 | S.missing.1)
    use.units   <- S %*% S.has.both * Eweights
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
    strat.sum.eweights <- t(S) %*% as.matrix(design@ElementWeights)
    strat.sum.eweights <- as.matrix(strat.sum.eweights)
    ## Horwitz Thompson-type assignment weights
    tx.wt <- ifelse(txclus.by.strat, nclus.by.strat/txclus.by.strat, 0)
    ctl.wt <-ifelse(ctlclus.by.strat, nclus.by.strat/ctlclus.by.strat, 0)
    tx.wt <- tx.wt/strat.sum.eweights   # appropriate HT weights for estimating, by stratum,
    ctl.wt <- ctl.wt/strat.sum.eweights # (total eweighted measurements)/(total eweights)
    ## factor in stratum weights
    tx.wt <- tx.wt * Swts
    ctl.wt <- ctl.wt * Swts
    ## now expand these up to match dimensions of data
    tx.wts <- as.vector(as.matrix(S %*% tx.wt))
    ctl.wts <- as.vector(as.matrix(S %*% ctl.wt))
    
    # ok, now that preliminaries are out of the way, compute some useful stuff.
    ## ratio estimates of means for a "domain" equal to intersection of
    ## treatment group with elements for which which the covariate is non-missing
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
##' totals up all the covariates
##' @param design DesignOptions
##' @return Aggregates
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

  element.weights <- as.matrix(t(C) %*% as.matrix(design@ElementWeights))
  dim(element.weights) <- NULL
  
  Eweights.tall <- design@ElementWeights * design@NotMissing
  Covariates <- as.matrix(t(C) %*% ifelse(Eweights.tall[,design@NM.Covariates, drop=FALSE],
                                          design@Covariates *
                                              Eweights.tall[,design@NM.Covariates, drop=FALSE],
                                          0)
                          )
  Eweights <- as.matrix(t(C) %*% Eweights.tall)
  Covariates <- ifelse(Eweights[,design@NM.Covariates, drop=FALSE] > 0,
                       Covariates/Eweights[,design@NM.Covariates, drop=FALSE],
                       0)
  NotMissing <- ifelse(matrix(element.weights>0, nrow(Eweights), ncol(Eweights)),
                       Eweights/element.weights, 0)
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
      ElementWeights = element.weights,
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
#' assumed to be nested within the stratifying variable. Extends DesignMatrix.
#'
#' Just as the covariates are presumed to have been aligned, the NotMissing weights carried
#' in this class are presumed to have been normalized: in each stratum, either all
#' weights are 0, or they've been rescaled to sum to 1. In order to handle different
#' NA patterns for different covariates, this normalization is done separately for each variable.
#'
#' In addition, unlike DesignOptions this class has no ElementWeights slot.
#' So the burden of representing differences in stratum
#' sizes falls entirely on the StrataWeightRatio slot. This slot has an entry for each unit,
#' representing ratio of user provided or specified stratum weight to h_b, the harmonic
#' mean of n_{tb} and n_{cb}, the counts of treatment and control clusters in stratum b.
#' (So h_b in no way reflects element weights or cluster sizes.)
#' 
#' @slot Z Logical indicating treatment assignment
#' @slot StrataMatrix A sparse matrix with n rows and s columns, with 1 if the unit is in that stratification
#' @slot StrataFactor Factor indicating strata
#' @slot StrataWeightRatio For each unit, ratio of stratum weight to h_b; but see Details.
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @slot OriginalVariables Look up table associating Covariates cols to terms in the calling formula, as in DesignMatrix
#' @slot Covariates Numeric matrix, as in DesignMatrix, except here we presume columns to have been aligned (stratum-centered)
#' @slot NotMissing matrix of element weights, normalized by stratum to sum to 1. (As otherwise, NAs are 0s)
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
alignDesignsByStrata <- function(design, post.align.transform = NULL) {

  stopifnot(inherits(design, "StratumWeightedDesignOptions")) # defensive programming

  vars   <- colnames(design@Covariates)
  strata <- names(design@StrataMatrices)

  Ewts  <- design@ElementWeights * design@NotMissing
  Covs <- ifelse(design@NotMissing[, design@NM.Covariates, drop = FALSE],
                 design@Covariates[, , drop = FALSE], 0)
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  Covs <- cbind(Covs, design@NotMissing[,NMcolperm])
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
    covars <- Covs[keep,,drop=FALSE]
    
    stopifnot(nlevels(ss)==1 ||
                  all(levels(ss)  %in% names(design@Sweights[[s]]$wtratio) ) )
    wtr.short <- design@Sweights[[s]]$wtratio[match(levels(ss),
                                                    names(design@Sweights[[s]]$wtratio),
                                                    nomatch=1L) # <-- this is to handle 
                                              ]                 # the unstratified case
    wtr <- wtr.short[ as.integer(ss) ]
    dim(wtr) <- NULL

    # align weighted observations within stratum by subtracting weighted stratum means
    covars.Sctr <- covars
    for (jj in 1L:ncol(covars)) 
        covars.Sctr[,jj] <-
            suppressWarnings( #throws singularity warning if covar is linear in S
                slm.wfit.csr( # see note in ./utils.R on why we use our own 
                    S, covars[,jj],               #`slm.wfit.csr` instead of `SparseM::slm.wfit`
                    weights=ewts[, covars.nmcols[jj], drop = TRUE])$residuals
                )

    if (!is.null(post.align.transform)) {
      # Transform the columns of covars.Sctr using the function in post.align.trans
      covars.Sctr.new <- apply(covars.Sctr, 2, post.align.transform)

      # Ensure that post.align.trans wasn't something that changes the size of covars.Sctr (e.g. mean).
      # It would crash later anyway, but this is more informative
      if (is.null(dim(covars.Sctr.new)) || !all(dim(covars.Sctr) == dim(covars.Sctr.new))) {
        stop("Invalid post.alignment.transform given")
      }
      ## The post alignment transform may have disrupted the stratum alignment.  So, recenter on stratum means
    for (jj in 1L:ncol(covars)) {
        covars.Sctr[,jj] <- suppressWarnings(
            slm.wfit.csr(S, covars.Sctr.new[,jj],
                         weights=ewts[, covars.nmcols[jj], drop = TRUE])$residuals
            )
    }
  }
    colnames(covars.Sctr) <- vars
    ewts.mns <- slm.fit.csr.fixed(S, ewts)$fitted
    ewts.mns <-  ifelse(ewts.mns==0, 1, ewts.mns)
    ewts.normed <- ewts/ewts.mns
      
      new("AlignedCovs",
          Z=as.logical(design@Z[keep]),
          StrataMatrix=S,
          StrataFactor=ss,
          StrataWeightRatio = wtr, #as extracted from the design
          Cluster           = factor(design@Cluster[keep]),
          OriginalVariables = origvars,
          Covariates        = covars.Sctr,
          NotMissing          = ewts.normed,
          TermLabels=design@TermLabels,
          Contrasts=design@Contrasts,
          NM.Covariates=covars.nmcols,
          NM.terms=design@NM.terms)
}
  lapply(strata, f)
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
    Eweights <- alignedcovs@NotMissing[,alignedcovs@NM.Covariates]
    Covs <- alignedcovs@Covariates
    
    n <- t(S) %*% S
    n.inv <- 1 / n
    n1 <- t(S) %*% zz
    n0 <- t(S) %*% (1 - zz)

    # dv is sample variance of treatment by stratum
    # set up 1/(n-1)
    tmp <- n
    tmp@ra <- 1 / (tmp@ra - 1)
    dv <- sparseToVec(S %*% tmp %*% (n1 - n.inv %*% n1^2)) * wtr^2
    
    tmat <- Covs * Eweights * #the sum statistic we're about to compute corresponds 
        wtr # to averaging within-stratum difference w/ stratum weights proportion to
    ## harmonic means of n_ts and n_cs.  To override w/user-designated weights, we
    ## factor in wtr, the "weight ratio" as previously reconstructed.
    ZtH <- S %*% n.inv %*% n1
    ssn <- sparseToVec(t(matrix(zz, ncol = 1) - ZtH) %*% tmat, column = FALSE)
    ssvar <- colSums(dv * tmat^2)
  
    ##  wtsum is the sum across strata of twice the harmonic mean of n1, n0 - we should rename it
    wtsum <- sum((n.inv %*% (n1 * n0))@ra) # (the ra slot is where SparseM keeps the non-zero values)

    zstat <- ifelse(ssvar <= .Machine$double.eps, NA_real_, ssn/sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

  scaled.tmat <- as.matrix(tmat * sqrt(dv))

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
  tcov <- crossprod(sqrt(dv) * tmat * (1 / wtsum))

  list(z = zstat, p = p, csq = csq , DF = DF, tcov = tcov)
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

