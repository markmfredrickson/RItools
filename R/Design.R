################################################################################
# Design Objects: communicate cluster, strata, treatment, and covariate information
################################################################################

#' Design S4 class
#'
#' @slot Z Logical indicating treatment assignment
#' @slot StrataMatrices This is a list of sparse matrices, each with n rows and s columns, with 1 if the unit is in that stratification
#' @slot StrataFrame Factors indicating strata (not the sparse matrices, as we use them in the weighting function)
#' @slot Cluster Factor indicating who's in the same cluster with who
#' @slot OriginalVariables Look up table to remind us which Covariates columns correspond to which provide variables
#' @slot Covariates Numeric matrix encoding variable values, analogous to a design matrix
#' @slot Eweights Matrix of element weights or cluster sums of them, with missing values contributing 0s.
setClass("Design",
         representation = list(
           Z                 = "logical",
           StrataMatrices    = "list", 
           StrataFrame       = "data.frame",  
           Cluster           = "factor",
           OriginalVariables = "character",
           Covariates        = "matrix",
           Eweights        = "matrix") 
         )

#
#
##' Create a Design object from a formula and some data
##'
##' The formula must have a left hand side that can be converted to a logical.
##'
##' On the RHS:
##'  - It may have at most one cluster() argument.
##'  - It may have one or more strata() arguments.
##'  - All other variables are considered covariates.
##'
##' NAs in a cluster() or strata() variable will be dropped.  NAs in covariates will be imputed.
##' 
##' @param fmla Formula
##' @param data Data
##' @param imputefn Function to impute
##' @param na.rm Should NA's removed?
##' @param include.NA.flags Flag NA's?
##' @return Design
makeDesign <- function(fmla, data, imputefn = median, na.rm = FALSE, include.NA.flags = TRUE) {

    eweights <- as.vector(model.weights(data))
    if (is.null(eweights))
        stop("makeDesign() expects its data arg to be a model frame containing weights")
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
    if (!(all(isGood))) {
      stop("In ", s, ", the following stratum levels did not include both treated and control units: ",
           paste(collapse = ", ", colnames(tbl)[!isGood]))
    }
  }

  ## OK! data looks good. Let's proceed to make the design object with covariate data

  data.fmla <- update(ts, paste("~", paste0(collapse = " - ", c(".", "1", vnames[str.idx]))))
  data.data <- model.frame(data.fmla, data, na.action = na.pass) #

  # knock out any levels that are not used
  fcts <- colnames(data.data)[sapply(data.data, is.factor)]
  for (f in fcts) {
    data.data[, f] <- factor(data.data[, f])
  }

  if (!na.rm) {
    # impute, possibly adding flags.
    data.data.imp <- naImpute(data.fmla, data.data, imputefn, include.NA.flags = include.NA.flags)
  } else {
    # who's missing entries in data.data
    idx <- !apply(data.data, 1, function(i) { any(is.na(i)) })
    data.data <- data.data[idx, ]
    data.data.imp <- data.data
    str.data <- str.data[idx, ]
  }

  # we want our own contrast function for the factors that expands each level to its own dummy
  tlbl <- names(data.data.imp)
  names(tlbl) <- as.character(tlbl)
  clist <- lapply(data.data.imp, function(x) {
    if (is.factor(x)) {
      structure(diag(nlevels(x)), dimnames = list(levels(x), levels(x)))
    } else {
      NULL
    }
  })
  clist <- clist[!sapply(clist, is.null)]

  data.mm         <- model.matrix(terms(data.data.imp), data.data.imp, contrasts.arg = clist)
  data.notmissing <- ifelse(is.na(model.matrix(terms(data.data), data.data, constrasts.arg = clist)),
                            0, eweights)

  # now we need to find if we added any NA flags
  toAdd <- dim(data.mm)[2] - dim(data.notmissing)[2]
  if (toAdd > 0) {
    xtra <- matrix(eweights, ncol = toAdd, nrow = dim(data.mm)[1])
    data.notmissing <- cbind(data.notmissing, xtra)
  }

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

  # create a look up table linking the model.matrix variables with the original variables

  originals <- attr(terms(data.fmla, data = data.data.imp), "term.labels")[attr(data.mm, "assign")]

  return(new("Design",
             Z                 = as.logical(as.numeric(Z) - 1),
             StrataMatrices    = strata.mats,
             StrataFrame       = strata.frame,
             Cluster           = factor(Cluster),
             Covariates        = data.mm,
             Eweights        = data.notmissing,
             OriginalVariables = originals))
}

# Add stratum weights to a design
#' Stratum Weighted Design
#'
#' @slot Sweights stratum weights
setClass("StratumWeightedDesign",
         representation = list(Sweights = "list"),
         contains = "Design")

## NB: Maybe later there'll a ClusterWeightedDesign class also,
## with methods for passing between the two.  ClusterWeighted would
## associate weights with each cluster, representing products of stratum
## weights and within-stratum HT-type weights.  Designs other than complete
## random assignment within blocks could be represented with such a structure.

##' Create stratum weights to be associated with a Design
##'
##' Apply weighting function to a Design by stratum, returning results in a format
##' suitable to be associated with an existing design and used in further calcs.
##' @param design Design
##' @param stratum.weights Stratum weights
##' @return data frame with a row for each stratum and two cols,
##' \code{sweights}  and \code{wtratio}. 
DesignWeights <- function(design, stratum.weights = harmonic) {
  stopifnot(inherits(design, "Design"))

  n.strata <- dim(design@StrataFrame)[2]
  strata.names <- colnames(design@StrataFrame)

  if (is.function(stratum.weights)) {
    swt.ls <- rep(list(stratum.weights), n.strata)
    names(swt.ls) <- strata.names
  }

  if (is.list(stratum.weights) & !all(strata.names %in% names(stratum.weights)))
    stop("list stratum.weights must have entry names matching those of stratifying factors")

  if (!is.list(stratum.weights) & !is.function(stratum.weights) & n.strata > 1)
    stop("stratum weights must be specified for each stratifying factor")

  if (!is.list(stratum.weights) & !is.function(stratum.weights)) {
    swt.ls <- list(stratum.weights)
    names(swt.ls) <- strata.names
  }

  if (is.list(stratum.weights) & !is.function(stratum.weights)) {
    swt.ls <- stratum.weights
    names(swt.ls) <- strata.names
  }


  ### change names here!

  wtlist <- list()
  for (nn in names(swt.ls)) {

    if (is.function(swt.ls[[nn]])) {
      sweights <-
        do.call(swt.ls[[nn]],
                args = list(data =
                    data.frame(Tx.grp = design@Z,
                               stratum.code = factor(design@StrataFrame[[nn]]),
                               design@Covariates,
                               check.names = FALSE)),
                envir=parent.frame())
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
                                  design@Covariates,
                                  check.names = FALSE))
    }
    hwts <- hwts/sum(hwts, na.rm=TRUE)

    wtratio <- unsplit(sweights/hwts, design@StrataFrame[[nn]], drop=TRUE)
    wtratio[is.na(wtratio)] <- 0
    wtlist[[nn]] <- list(sweights=sweights,wtratio=wtratio)
    NULL
  }
  wtlist
}

##' Generate Descriptives
##'
##' Use a design object to generate descriptive statistics that ignore
##' clustering. Stratum weights are respected.
##' @param design A Design
##' @param covariate.scaling Scale estimates for covs, to use instead of internally calculated pooled SDs
##' @return Descriptives
designToDescriptives <- function(design, covariate.scaling = NULL) {
  stopifnot(inherits(design, "Design")) # defensive programming
  if (!is.null(covariate.scaling)) warning("Non-null 'covariate.scaling' currently being ignored")
  vars   <- colnames(design@Covariates)
  stratifications <- colnames(design@StrataFrame)

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
    if (inherits(design, "StratumWeightedDesign") & length(stratlevs)>1)
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

    S.missing.0 <- as.matrix((t(ZZ) %*% design@Eweights)) == 0
    S.missing.1 <- as.matrix((t(WW) %*% design@Eweights)) == 0
    S.has.both  <- !(S.missing.0 | S.missing.1)
    use.units   <- S %*% S.has.both * design@Eweights

    X.use  <- design@Covariates * use.units
    X2.use <- design@Covariates^2 * use.units

    n1 <- t(use.units) %*% Z
    n0 <- t(use.units) %*% (1 - Z)
    
    ## Now calculate assignment/stratum weights
    cluster.representatives <- !duplicated(design@Cluster)
    txcts <- table(as.logical(Z[cluster.representatives]),
                   design@StrataFrame[cluster.representatives, s])
    stratcts <- txcts["TRUE",,drop=FALSE] + txcts["FALSE",,drop=FALSE]
    ## HT-type assignment weights
    tx.wt <- ifelse(txcts["TRUE",,drop=FALSE],stratcts/txcts["TRUE",,drop=FALSE], 0)
    ctl.wt <-ifelse(txcts["FALSE",,drop=FALSE],stratcts/txcts["FALSE",,drop=FALSE], 0)
    ## multiplying through for assignment/stratum weights
    tx.wt <- t(tx.wt[,stratlevs,drop=FALSE]) * Swts
    ctl.wt <- t(ctl.wt[,stratlevs, drop=FALSE]) * Swts
    ## now expand these up to match dimensions of data
    tx.wts <- as.vector(as.matrix(S %*% tx.wt))
    ctl.wts <- as.vector(as.matrix(S %*% ctl.wt))
    
    # ok, now that preliminaries are out of the way, compute some useful stuff.
    wtsum.tx <-  t(use.units * tx.wts) %*% Z
    treated.avg <- t(X.use *tx.wts) %*% Z / wtsum.tx

    wtsum.ctl <- t(use.units * ctl.wts) %*% (1 - Z)
    control.avg <- t(X.use * ctl.wts) %*% (1 - Z) / wtsum.ctl

    var.1 <- (t(X2.use *tx.wts) %*% Z - wtsum.tx * treated.avg^2) / wtsum.tx
    var.1 <- var.1 * n1/(n1 - 1)
    var.0 <- (t(X2.use * ctl.wts) %*% (1 - Z) - wtsum.ctl * control.avg^2) / wtsum.ctl
    var.0 <- var.0* n0/(n0 - 1)

    pooled <- sqrt((var.1 + var.0) / 2)

    adjustedDifference    <- treated.avg - control.avg
    standardizedDifference <- adjustedDifference / pooled

    ans[, , s] <- c(control.avg, treated.avg, standardizedDifference, adjustedDifference, pooled)
  }

  return(ans)
}

##' Aggregate Design
##'
##' totals up all the covariates
##' @param design Design
##' @return Aggregates
aggregateDesign <- function(design) {
  n.clusters <- nlevels(design@Cluster)

  if (n.clusters == length(design@Cluster)) {
    return(design)
  }

  dupes <- duplicated(design@Cluster)
  Z <- design@Z[!dupes]
  StrataFrame <- design@StrataFrame[!dupes,, drop = FALSE]
  Cluster <- design@Cluster[!dupes]

  C <- SparseMMFromFactor(design@Cluster)
  
  Eweights <- as.matrix(t(C) %*% design@Eweights)
  Covariates <- as.matrix(t(C) %*% (design@Covariates * design@Eweights))
  Covariates <- ifelse(Eweights>0, Covariates/Eweights, 0)
  
  StrataMatrices <- lapply(design@StrataMatrices, function(S) {
    tmp <- t(C) %*% S
    tmp@ra <- rep(1, length(tmp@ra))
    return(tmp)
  })

  # colnames(Covariates) <- c("cluster.size", colnames(design@Covariates))
  colnames(Covariates)   <- colnames(design@Covariates)

  new("Design",
      Z = Z,
      StrataMatrices = StrataMatrices,
      StrataFrame = StrataFrame,
      Cluster = Cluster,
      Eweights = Eweights,
      Covariates = as.matrix(Covariates),
      # OriginalVariables = c("Cluster Size", design@OriginalVariables))
      OriginalVariables = design@OriginalVariables)

}

##' Align Design by Strata
##'
##' TODO: Can we minimize the amount of weighted stuff that needs to get pushed through?
##' Eg. when creating the weighted design, we mulitply all covariates by the weighting scheme
##' right away, so we don't need to track it explicitly later.
##' observe all the multiplication of tmat by swt$wtradio throughout. This seems redundant
##' @param design Design
##' @param post.align.transform A post-align transform
##' @return List
alignDesignByStrata <- function(design, post.align.transform = NULL) {

  stopifnot(inherits(design, "StratumWeightedDesign")) # defensive programming

  vars   <- colnames(design@Covariates)
  strata <- names(design@StrataMatrices)

  # we can't return an array because different stratifications will have varying numbers
  # of strata levels. A list is more flexible here, but less structured.
  ans <- list()

  for (s in strata) {
    ss <- design@StrataFrame[, s]
    keep <- !is.na(ss)
    ss <- ss[keep]
    S <- SparseMMFromFactor(ss)
    Z <- as.numeric(design@Z[keep])
    n <- t(S) %*% S
    n.inv <- 1 / n
    n1 <- t(S) %*% Z
    n0 <- t(S) %*% (1 - Z)

    Covs <- design@Covariates[keep, , drop = FALSE]
    Ewts <- design@Eweights[keep, , drop = FALSE]
    wtr <- design@Sweights[[s]]$wtratio[keep]

    # see note in ./utils.R on why we use this function instead of SparseM's version
    # align weighted observations within stratum by subtracting weighted stratum means
    Covs.Sctr <- Covs
    for (jj in 1:ncol(Covs))
      Covs.Sctr[,jj] <- slm.wfit.csr(S, Covs[,jj], weights=Ewts[,jj])$residuals

    # dv is sample variance of treatment by stratum
    # set up 1/(n-1)
    tmp <- n
    tmp@ra <- 1 / (tmp@ra - 1)

    dv <- sparseToVec(S %*% tmp %*% (n1 - n.inv %*% n1^2)) * wtr^2


    if (!is.null(post.align.transform)) {
      # Transform the columns of Covs.Sctr using the function in post.align.trans
      Covs.Sctr.new <- apply(Covs.Sctr, 2, post.align.transform)

      # Ensure that post.align.trans wasn't something that changes the size of Covs.Sctr (e.g. mean).
      # It would crash later anyway, but this is more informative
      if (is.null(dim(Covs.Sctr.new)) || !all(dim(Covs.Sctr) == dim(Covs.Sctr.new))) {
        stop("Invalid post.alignment.transform given")
      }
      ## The post alignment transform may have disrupted the stratum alignment.  So, recenter on stratum means
    for (jj in 1:ncol(Covs)) {
      Covs.Sctr[,jj] <- slm.wfit.csr(S, Covs.Sctr.new[,jj], weights=Ewts[,jj])$residuals }

  }

    
    tmat <- Covs.Sctr * Ewts * wtr
    ZtH <- S %*% n.inv %*% n1
    ssn <- sparseToVec(t(matrix(Z, ncol = 1) - ZtH) %*% tmat, column = FALSE)
    ssvar <- colSums(dv * tmat^2)
  
    ##  wtsum is the sum across strata of twice the harmonic mean of n1, n0 - we should rename it
    wtsum <- sum((n.inv %*% (n1 * n0))@ra) # the ra slot is where SparseM keeps the non-zero values)


    # save everything as we drop some of the observations and we need all the dims/etc to line up
    ans[[s]] <- list(zz = Z,
                     tmat = tmat, # we can make this dense, chances are the zeros are not especially common
                     ssn = ssn,
                     ssvar = ssvar,
                     dv = dv,
                     wtsum = wtsum)
  }
  return(ans)
}

# I'd prefer to have a better API here, but right now, just trying to get compatability with old xBalance set up
# e.g. something that is a list of strata with a given structure, rather than just a list.
##' Align to Inferentials
##'
##' @param zz zz
##' @param tmat tmat
##' @param ssn ssn
##' @param ssvar ssvar
##' @param dv dv
##' @param wtsum wtsum
##' @return list
alignedToInferentials <- function(zz, tmat, ssn, ssvar, dv, wtsum) {
  z <- ifelse(ssvar <= .Machine$double.eps, 0, ssn/sqrt(ssvar))
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  scaled.tmat <- as.matrix(tmat * sqrt(dv))

  pst.svd <- try(svd(scaled.tmat, nu=0))

  if (inherits(pst.svd, 'try-error')) {
    pst.svd <- propack.svd(scaled.tmat)
  }

  Positive <- pst.svd$d > max(sqrt(.Machine$double.eps) * pst.svd$d[1], 0)
  Positive[is.na(Positive)] <- FALSE # JB Note: Can we imagine a situation in which we dont want to do this?

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

  list(z = z, p = p, csq = csq , DF = DF, tcov = tcov)
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

