#' @include ModelMatrixPlus.R
NULL

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

  ## convert any character columns to factors
  for (cs in colnames(data.data)[sapply(data.data, is.character)]) {
    data.data[, cs] <- as.factor(data.data[, cs])
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
  colnames(tmp) <- gsub(colnames(tmp), pattern = "survival::strata\\((.*)\\)", replacement = "\\1")
  strata.frame <- data.frame(lapply(tmp, factor), check.names = FALSE)

  return(new("DesignOptions", desm,
             Z           = as.logical(as.numeric(Z) - 1), #b/c it was built as factor
             StrataFrame = strata.frame,
             Cluster     = factor(Cluster)))
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
  covars <- ifelse(is.na(design), 0, design)

  ## Tack NM cols onto covars, but with intercept col listed last
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  covars <- cbind(covars, design@NotMissing[, NMcolperm])
  vars <- c(colnames(design),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
  colnames(covars)   <-  vars
  covars.nmcols <- c(pmax(1L, design@NM.Covariates), rep(1L, k.NM ) )

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

  Uweights <- design@UnitWeights * design@NotMissing

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

    ## strat.sum.uweights *(txclus.by.strat / nclus.by.strat)
    tx.wt <- ifelse(txclus.by.strat, 1, 0)
    ctl.wt <-ifelse(ctlclus.by.strat, txclus.by.strat/ctlclus.by.strat, 0)

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
    var.1 <- (t(X2.use *tx.wts) %*% Z - wtsum.tx * treated.avg^2) / wtsum.tx
    var.1 <- var.1 * ifelse(n1>1, n1/(n1 - 1), 0)
    var.0 <- (t(X2.use * ctl.wts) %*% (1 - Z) - wtsum.ctl * control.avg^2) / wtsum.ctl
    var.0 <- var.0* ifelse(n0>1, n0/(n0 - 1), 0)

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
  # To align w/ this `C`, everything to be returned 
  # needs to align w/ `levels(Cluster)`, not with 
  # `Cluster` itself. So,
  Z <- Z[levels(Cluster)]
  Cluster  <- as.factor(levels(Cluster))
  StrataFrame  <-
      StrataFrame[match(levels(Cluster), as.character(Cluster)),
                  , drop=FALSE]
  
  unit.weights <- as.matrix(C_transp %*% as.matrix(design@UnitWeights))

  dim(unit.weights) <- NULL
  names(unit.weights) <- levels(Cluster)

  Uweights.tall <- design@UnitWeights * design@NotMissing
  Covariates <- as.matrix(C_transp %*% ifelse(Uweights.tall[,pmax(1L,design@NM.Covariates), drop=FALSE],
                                          design *
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
    colnames(Covariates)   <- colnames(design)
    row.names(Covariates) <- levels(Cluster)

  new("DesignOptions", Covariates,
      Z = Z,
      StrataFrame = StrataFrame,
      Cluster = Cluster,
      UnitWeights = unit.weights,
      NotMissing = NotMissing,
      OriginalVariables = design@OriginalVariables,
      TermLabels=design@TermLabels,
      Contrasts=design@Contrasts,
      NM.Covariates=design@NM.Covariates,
      NM.terms=design@NM.terms)
}



### Putting these two classes here to quiet warning about missing class.
### A better solution will use collation to get the right file order.

## A set of virtual classes to be used with specific kinds of designs
setClass("RandomizedDesign")

## Stratified design of x units into k strata with fixed numbers of units and
## treated per strata
##
## @slot Unit A n by k (sparse) matrix with a single 1 in each row. Rows
##   indicated the randomized units (clusters in cluster randomized trials).
##   Columns indicate stratum membership.
## @slot Counts A k vector of units in each stratum
## @slot Treated A k vector of the number treated in each stratum.
## @slot Weights A k vector of strata weights.
## @slot JCoefs A n vector of precomputes coefficients to produce J = Z (1 - P(Z = 1)) / (P(Z = 1) P(Z = 0))
setClass("StratifiedDesign", contains = "RandomizedDesign",
         slots = c(
             Units = "matrix.csr",
             Count = "integer",
             Treated = "integer",
             Weights = "numeric",
             JCoefs = "numeric"
         ))

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
#' specified stratum weight to the product of h_b (the harmonic mean of n_{tb} and
#' n_{cb}, the counts of treatment and control clusters in stratum b) with bar-w_b,
#' (the arithmetic mean of aggregated cluster weights within that stratum). It can also
#' be the numeric vector 1, without names, meaning the intended weight ratio is always 1.
#' 
#' @slot Covariates Numeric matrix, as in ModelMatrixPlus, except: will include NM columns; all columns presumed to have been stratum-centered (aligned)
#' @slot Z Logical indicating treatment assignment
#' @slot Design A StratafiedDesign object that captures the stratification and weighting 
#' @keywords internal
setClass("CovsAlignedToADesign",
         slots                 =
             c(Covariates      = "matrix",
             Z                 = "logical",
             Design            = "StratifiedDesign"
             )
         )



## Useful methods for RandomizedDesigns to implement:
## toJ: given a Z, return (Z - P(Z = 1)) / (P(Z = 1)P(Z = 0))
## squared_diff_covariance, return the Cov(T^2) k by k matrix

## Prepare the J vector from M = J' W J
##
## @param z A treatment assignment vector (length n).
## @return A vector of length n
toJ <- function(x, z) { UseMethod("toJ") }

## Return the Cov(T^2) matrix of rotated variable totals T_k
##
## We can write the Mahalanobis distance as \sum_{i = 1}^k T_i^2 where
## T = (xV)'J
## where xV is the design rotated covariates
##
## The variance of the Mahalanobis distance is then 1' Cov(T^2) 1
t_squared_covariance <- function(design, covariates) {
    UseMethod("t_squared_covariance")
}

## Return the Cov(S^2) matrix of variable totals S_k
euclidean_squared_covariance <- function(design, covariates) {
    UseMethod("euclidean_squared_covariance")
}

setClass("IndependentRandomizationDesign", contains = "RandomizedDesign",
         slots = c(InclusionProbabilities = "numeric"))

## Combines a design with a covariate matrix, modifying the covariate matrix by multiplying it 
setClass("DesignRotatedCovariates", contains = "matrix",
         slots = c(Design = "RandomizedDesign", Rotation = "matrix"))



## Modify covariates to have be multiplied by the inverse of the square root inverse covariance matrix of the mahalanobis statistic
##
## @param design A RandomizedDesign object
## @param x A covariate matrix
## @return a DesignRotatedCovariates object (subclass of matrix)
rotate_covariates <- function(design, x) {
    UseMethod("rotate_covariates")
}

## Return T = (xV)'J starting from a treatment assignment vector.
toT <- function(rotated, z) {
    stop("not implemented yet")
}



# apply this & pass through en route to svd
#' @method scale DesignOptions
scale.DesignOptions  <- function(x, center=TRUE, scale=TRUE)
{
    stopifnot(is(x, "DesignOptions"))
    refstrat  <- which(colnames(x@StrataFrame)=="--")
    refstrat  <- if (length(refstrat)==0) 1L else refstrat[1]
    x@StrataFrame <- x@StrataFrame[refstrat]
    wtsum  <- sum(x@UnitWeights[complete.cases(x@StrataFrame)])
    trans  <- if (is(center, "function")) center else NULL
    aligned  <- alignDesignsByStrata(x, post.align.transform=trans)
    aligned_covs  <- aligned[[1]]
    if (scale)
        {
    scales  <- .colSums(aligned_covs@Covariates^2,
                        nrow(aligned_covs@Covariates), ncol(aligned_covs@Covariates))
    scales  <- scales/wtsum
    scales[scales<.Machine$double.eps^.5]  <- 1
    sweep(aligned_covs@Covariates, 2L, scales, "/", check.margin=FALSE)
        } else aligned_covs@Covariates
    }

##' Align DesignOptions by Strata
##'
##' @param design DesignOptions
##' @param post.align.transform A post-align transform
##' @return list List of `CovsAlignedToADesign` objects
##' @keywords internal
##'
alignDesignsByStrata <- function(design, post.align.transform = NULL) {


  vars   <- colnames(design)
  stratifications <- colnames(design@StrataFrame)

  Ewts  <- design@UnitWeights * design@NotMissing
  Covs <- ifelse(design@NotMissing[, pmax(1L,design@NM.Covariates), drop = FALSE],
                 design[, , drop = FALSE], 0)
  k.Covs <- ncol(Covs)
  NMcolperm <- if ( (k.NM <- ncol(design@NotMissing)) >1) c(2L:k.NM, 1L) else 1
  Covs <- cbind(Covs, 0+design@NotMissing[,NMcolperm])
  vars  <- c(colnames(design),  paste0("(", colnames(design@NotMissing)[NMcolperm], ")") )
  covars.nmcols <- c(design@NM.Covariates, rep(1L, k.NM ) )
  origvars <- match(colnames(design@NotMissing), design@TermLabels, nomatch=0L)
  origvars <- c(design@OriginalVariables, origvars)

  # we can't return an array because different stratifications will have varying numbers
  # of strata levels. A list is more flexible here, but less structured.
  f <- function(s){
    ss <- design@StrataFrame[, s]
    keep <- !is.na(ss)
    ss <- ss[keep]
    zz <- design@Z[keep]


    ewts <- Ewts[keep,,drop=FALSE]
    non_null_record_wts <- ewts[,1L,drop=TRUE]
    NM <- design@NotMissing[keep,,drop=FALSE]
    NMCovs <- design@NM.Covariates
    covars <- Covs[keep,,drop=FALSE]

    stratified <- create_stratified_design(ss, z = zz)

    #### End move

    # align weighted observations within stratum by subtracting weighted stratum means
    covars.Sctr <- covars
    for (jj in 1L:ncol(covars))
    {
        covars.Sctr[,jj] <- if (all(covars[,jj]==covars[1L,jj])) 0 else
            suppressWarnings( #throws singularity warning if covar is linear in S
                slm.wfit.csr( # see note in ./utils.R on why we use our own
                    stratified@Units, covars[,jj],               #`slm.wfit.csr` instead of `SparseM::slm.wfit`
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
            slm.wfit.csr(stratified@Units, covars.Sctr.new[,1L:k.Covs, drop=FALSE],
                         weights=non_null_record_wts)$residuals
        )
    }
    colnames(covars.Sctr) <- vars
    covars.Sctr  <- covars.Sctr * non_null_record_wts

      new("CovsAlignedToADesign",
          Covariates = covars.Sctr,
          Z = zz,
          Design = stratified
          )
}
  sapply(stratifications, f, simplify = FALSE, USE.NAMES = TRUE)
}

##' @title Adjusted & combined differences as in Hansen & Bowers (2008)
##' @param alignedcovs A CovsAlignedToADesign object
##' @return list with components:
##' \describe{
##'   \item{z}{First item}
##'   \item{p}{Second item}
##'   \item{Msq}{Squared Mahalanobis distance of combined differences from origin}
##'   \item{DF}{degrees of freedom}
##'   \item{adj.mean.diffs}{Vector of sum statistics z'x-tilde, where x-tilde is the unit- and stratum-weighted covariate, with stratum centering.  This differs from the adjusted difference vector of Hansen & Bowers (2008) by a constant of proportionality.}
##'   \item{tcov}{Matrix of null covariances of Z'x-tilde vector, as above.}
##' }
##' @references Hansen, B.B. and Bowers, J. (2008), ``Covariate
##'   Balance in Simple, Stratified and Clustered Comparative
##'   Studies,'' \emph{Statistical Science} \bold{23}.
##' @seealso \code{\link{balanceTest}}, \code{\link{alignDesignsByStrata}}
##' @importMethodsFrom SparseM diag
##' @keywords internal
HB08 <- function(alignedcovs) {

    ## appropriately weighted sum of each requested variable
    ssn <- manifest_variable_sums(alignedcovs@Design, alignedcovs@Covariates, alignedcovs@Z)

    tcov <- manifest_variable_covariance(alignedcovs@Design, alignedcovs@Covariates)

    ssvar <- diag(tcov)

    zero_variance  <- (ssvar <= .Machine$double.eps)
    zstat <- ifelse(zero_variance, NA_real_, ssn/sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)


    cov_minus_.5 <-
        XtX_pseudoinv_sqrt(mat=tcov[!zero_variance, !zero_variance, drop=FALSE],
                           mat.is.XtX = TRUE)
    mvz <- drop(crossprod(ssn[!zero_variance], cov_minus_.5))
    csq <- drop(crossprod(mvz))
    DF <- ncol(cov_minus_.5)

    list(z = zstat, p = p, Msq = csq , DF = DF,
         adj.mean.diffs=ssn, tcov = tcov)
}

##' @title Hansen & Bowers (2008) inferentials 2016 [81e3ecf] version
##' @param alignedcovs A CovsAlignedToADesign object
##' @return list, as in \code{\link{HB08}}
##' @keywords internal
HB08_2016 <- function(alignedcovs) {
    zz <- as.numeric(alignedcovs@Z)
    S <- alignedcovs@Design@Units

    Covs <- alignedcovs

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

    n1_over_n <- S %*% n.inv %*% n1
    ssn <- sparseToVec(t(matrix(zz, ncol = 1) - n1_over_n) %*% Covs, column = FALSE)
    names(ssn) <- colnames(Covs)

    scaled.Covs <- as.matrix(Covs * sqrt(dv))
    tcov <- crossprod(scaled.Covs)
    ssvar <- diag(tcov)

    zero_variance  <- (ssvar <= .Machine$double.eps)
    zstat <- ifelse(zero_variance, NA_real_, ssn/sqrt(ssvar))
    p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

    ## moving forward, we'll do without those sum statistics that have 0 null variation.
    x_tilde <- Covs[,zero_variance, drop=FALSE]
    scaled.x_tilde <- scaled.Covs[,ssvar > .Machine$double.eps, drop=FALSE]
    ssn  <- ssn[ssvar > .Machine$double.eps]

    cov_minus_.5 <- XtX_pseudoinv_sqrt(scaled.x_tilde)
    mvz <- drop(crossprod(ssn, cov_minus_.5))
    csq <- drop(crossprod(mvz))
    DF <- ncol(cov_minus_.5)

    list(z = zstat, p = p, Msq = csq , DF = DF,
       adj.mean.diffs=ssn, tcov = tcov)
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
