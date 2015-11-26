################################################################################
# Design Objects: communicate cluster, strata, treatment, and covariate information
################################################################################

setClass("Design",
         representation = list(
           Z                 = "logical",
           Strata            = "data.frame",
           Cluster           = "factor",
           OriginalVariables = "character",
           Covariates        = "matrix",
           NotMissing        = "matrix") # 1 is the value in covariates was not missing, 0 if it was imputed
         )
             

# Create a Design object from a formula and some data
#
# The formula must have a left hand side that can be converted to a logical.
#
# On the RHS:
#  - It may have at most one cluster() argument.
#  - It may have one or more strata() arguments.
#  - All other variables are considered covariates.
#
# NB: should we make this more like glm() to pick up environment? Probably
makeDesign <- function(fmla, data, imputefn = median, na.rm = FALSE, include.NA.flags = TRUE) {
  ts <- terms(fmla, data = data, specials = c("cluster", "strata"))

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
  str.data <- model.frame(str.tms, data = data, na.action = na.pass)
  

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
  data.notmissing <- 1 - is.na(model.matrix(terms(data.data), data.data, constrasts.arg = clist))

  # now we need to find if we added any NA flags
  toAdd <- dim(data.mm)[2] - dim(data.notmissing)[2]
  if (toAdd > 0) {
    xtra <- matrix(1, ncol = toAdd, nrow = dim(data.mm)[1])
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
  Strata <- data.frame(lapply(tmp, factor), check.names = FALSE)
  
  # create a look up table linking the model.matrix variables with the original variables
  
  originals <- attr(terms(data.fmla, data = data.data.imp), "term.labels")[attr(data.mm, "assign")]

  return(new("Design",
             Z = as.logical(as.numeric(Z) - 1),
             Strata = Strata,
             Cluster = factor(Cluster),
             Covariates = data.mm,
             NotMissing = data.notmissing,
             OriginalVariables = originals))
}

# Add stratum weights to a design
setClass("WeightedDesign",
         representation = list(Weights = "list"),
         contains = "Design")

weightedDesign <- function(design, stratum.weights = harmonic, normalize.weights = TRUE) {
  stopifnot(inherits(design, "Design"))
  
  n.strata <- dim(design@Strata)[2]
  strata.names <- colnames(design@Strata)
  
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
                               stratum.code = factor(design@Strata[[nn]]),
                               design@Covariates,
                               check.names = FALSE)),
                envir=parent.frame())
    } else {
      if (!is.numeric(swt.ls[[nn]]))
        stop("stratum.weights must be an expression or numeric vector")

      if (is.null(names(swt.ls[[nn]])))
        stop ("if stratum.weights is a vector, must have names")

      if (!(all(levels(factor(design@Strata[[nn]])) %in% names(swt.ls[[nn]])) ))
        stop("if stratum.weights is a vector, must have a name for each stratum")

      sweights <- swt.ls[[nn]][levels(factor(design@Strata[[nn]]))]
    }

    if (all(is.na(sweights)))
      stop(paste("All stratum weights NA (strat.",nn,")."))
    if (any(is.na(sweights))) {
      sweights[is.na(sweights)] <- 0
      warning(paste("NAs in stratum weights (",nn," strat.); to be interpreted as 0s.", sep=""))
    }
    if (any(sweights<0))
      stop("stratum weights must be nonnegative")

    if (normalize.weights)
      sweights <- sweights/sum(sweights, na.rm=TRUE)

    if (identical(harmonic, swt.ls[[nn]])) {
      hwts <- sweights
    } else {
      hwts <- harmonic(data.frame(Tx.grp = design@Z,
                                  stratum.code=factor(design@Strata[[nn]]),
                                  design@Covariates,
                                  check.names = FALSE))
    }
    hwts <- hwts/sum(hwts, na.rm=TRUE)

    wtratio <- unsplit(sweights/hwts, design@Strata[[nn]], drop=TRUE)
    wtratio[is.na(wtratio)] <- 0
    wtlist[[nn]] <- list(sweights=sweights,wtratio=wtratio)
    NULL
  }
  wtlist

  design <- as(design, "WeightedDesign")
  design@Weights <- wtlist
  return(design)
}

# Use a design object to generate descriptive statistics that ignore clustering.
# Stratum weights are respected.
designToDescriptives <- function(design, covariate.scaling = TRUE) {
  stopifnot(inherits(design, "Design")) # defensive programming

  vars   <- colnames(design@Covariates)
  strata <- colnames(design@Strata)

  ans <- array(NA,
               dim = c(length(vars), 5, length(strata)),
               dimnames = list(
                   "vars" = vars,
                   "stat" = c("Control", "Treatment", "std.diff", "adj.diff", "pooled.sd"),
                   "strata" = strata))
  for (s in strata) {
    
    S <- SparseMMFromFactor(design@Strata[, s])
    Z <- as.numeric(design@Z)
    ZZ <- S * Z 
    WW <- S * (1 - Z)
    
    S.missing.0 <- as.matrix((t(ZZ) %*% design@NotMissing)) == 0 
    S.missing.1 <- as.matrix((t(WW) %*% design@NotMissing)) == 0 
    S.has.both  <- !(S.missing.0 | S.missing.1)
    use.units   <- S %*% S.has.both * design@NotMissing

    X.use  <- design@Covariates * use.units
    X2.use <- X.use^2 

    n1 <- t(use.units) %*% Z
    n0 <- t(use.units) %*% (1 - Z)

    ETT <- S %*% t(ZZ) %*% use.units

    # ok, now that preliminaries are out of the way, compute some useful stuff.
    treated.avg <- t(X.use) %*% Z / n1

    n0.ett <- t(ETT) %*% (1 - Z)
    control.avg <- t(X.use * ETT) %*% (1 - Z) / n0.ett

    var.1 <- (t(X2.use) %*% Z - n1 * treated.avg^2) / (n1 - 1)
    var.0 <- (t(X2.use * ETT) %*% (1 - Z) - n0.ett * control.avg^2) / (n0.ett - 1)

    pooled <- sqrt((var.1 + var.0) / 2)
    
    adjustedDifference    <- treated.avg - control.avg
    standardizedDifference <- adjustedDifference / pooled

    ans[, , s] <- c(control.avg, treated.avg, standardizedDifference, adjustedDifference, pooled)
  }

  return(ans)
}

## <description>
## Turn a factor variable into a sparse matrix of 0's and 1's, such that if observation i
## has the jth level then there is a 1 at position (i,j) (but nowhere else in row i).
## <details>
## NA's give rise to rows with no 1s.
## As the result is only meaningful in the context of the SparseM package,
## function requires that SparseM be loaded.
## @title Sparse matrix dummy coding of a factor variable (omitting the intercept)
## @param thefactor Factor variable, or object inheriting from class factor
## @return Sparse csr matrix the columns of which are dummy variables for levels of thefactor
## @import SparseM
## @export
## @author Ben Hansen
## @examples
## sparse_mod_matrix <-  SparseMMFromFactor(iris$Species)
## mod_matrix <- model.matrix(~Species-1, iris)
## all.equal(as.matrix(sparse_mod_matrix),
##           mod_matrix, check.attributes=FALSE)
SparseMMFromFactor <- function(thefactor) {
  stopifnot(inherits(thefactor, "factor"))
  theNA <- ##if (inherits(thefactor, "optmatch")) !matched(thefactor) else
    is.na(thefactor)

  if (all(theNA)) stop("No non-NA's in thefactor") else {
    if (any(theNA) && !inherits(thefactor, "optmatch")) warning("NA's found in thefactor.")
  }

  nlev <- nlevels(thefactor)
  nobs <- length(thefactor)
  theint <- as.integer(thefactor)
  if (any(theNA)) theint[theNA] <- 1L#nlev + 1L:sum(theNA)
  new("matrix.csr",
      ja=theint,
      ia=1L:(nobs+1L),
      ra=(1L-theNA),
      dimension = c(nobs, nlev) #+sum(theNA)
      )
}
# Aggregated Design totals up all the covariates 

aggregateDesign <- function(design) {
  n.clusters <- nlevels(design@Cluster)

  if (n.clusters == length(design@Cluster)) {
    return(design)
  }

  # it seemed like a good idea to include cluster size, but this leads to issues later with z-values and p-values for the descriptives
  # keeping the idea commented out to make it easier to resurrect later
  # Covariates <- matrix(NA, nrow = n.clusters, ncol = ncol(design@Covariates) + 1) # the extra col will be for cluster counts
  Covariates <- matrix(NA, nrow = n.clusters, ncol = ncol(design@Covariates)) 

  dupes <- duplicated(design@Cluster)
  Z <- design@Z[!dupes]
  Strata <- design@Strata[!dupes,, drop = FALSE]
  Cluster <- design@Cluster[!dupes]

  for (i in 1:length(Cluster))  {
    cname <- Cluster[i]
    subcovs <- design@Covariates[design@Cluster == cname, , drop = FALSE]

    # Covariates[i, ] <- c(dim(subcovs)[1], colSums(subcovs))
    Covariates[i, ] <- colSums(subcovs)
  }


  # colnames(Covariates) <- c("cluster.size", colnames(design@Covariates))
  colnames(Covariates)   <- colnames(design@Covariates)

  new("Design",
      Z = Z,
      Strata = Strata,
      Cluster = Cluster,
      Covariates = as.matrix(Covariates),
      # OriginalVariables = c("Cluster Size", design@OriginalVariables))
      OriginalVariables = design@OriginalVariables)
                       
}

# TODO: Can we minimize the amount of weighted stuff that needs to get pushed through?
# Eg. when creating the weighted design, we mulitply all covariates by the weighting scheme
# right away, so we don't need to track it explicitly later.
# observe all the multiplication of tmat by swt$wtradio throughout. This seems redundant
alignDesignByStrata <- function(design, post.align.transform = NULL) {

  stopifnot(inherits(design, "WeightedDesign")) # defensive programming

  vars   <- colnames(design@Covariates)
  strata <- colnames(design@Strata)

  # we can't return an array because different stratifications will have varying numbers
  # of strata levels. A list is more flexible here, but less structured.
  ans <- list()

  for (s in strata) {

    S <- SparseMMFromFactor(design@Strata[, s])
    Z <- as.numeric(design@Z)
    n <- t(S) %*% S
    n.inv <- 1 / n
    n1 <- t(S) %*% Z
    n0 <- t(S) %*% (1 - Z)

    ZtH <- S %*% n.inv %*% n1
    ssn <- sparseToVec(t(matrix(Z, ncol = 1) - ZtH) %*% (design@Covariates * design@Weights[[s]]$wtratio), column = FALSE)

    wtsum <- sum((n.inv %*% (n1 * n0))@ra) # the ra slot is where SparseM keeps the non-zero values)

    msmn <- as.matrix(S %*% n.inv %*% t(S) %*% design@Covariates) # this will be mostly dense, so just cast back to a standard matrix

    tmat <- (design@Covariates - msmn)

    # dv is sample variance of treatment by stratum
    # set up 1/(n-1)
    tmp <- n
    tmp@ra <- 1 / (tmp@ra - 1)

    dv <- sparseToVec(S %*% tmp %*% (n1 - n.inv %*% n1^2))
    ssvar <-  colSums(dv * design@Weights[[s]]$wtratio^2 * tmat^2)
    
    if (!is.null(post.align.transform)) {
      # Transform the columns of tmat using the function in post.align.trans
      tmat.new <- apply(tmat, 2, post.align.transform)

      # Ensure that post.align.trans wasn't something that changes the size of tmat (e.g. mean).
      # It would crash later anyway, but this is more informative
      if (is.null(dim(tmat.new)) || !all(dim(tmat) == dim(tmat.new))) {
        stop("Invalid post.alignment.transform given")
      }
      ## recenter on stratum means
      tmat <- tmat.new
      msmn <- as.matrix(S %*% n.inv %*% t(S) %*% tmat)
      tmat <- tmat - msmn
      tmat <- tmat * design@Weights[[s]]$wtratio

      # Recalculate on transformed data the difference of treated sums and their null expectations
      # (NB: since tmat has just been recentered,
      ssn <- t(Z) %*% tmat
      ssvar <- colSums(dv * tmat^2)
    } else {
      tmat <- tmat * design@Weights[[s]]$wtratio
    }

    # save everything as we drop some of the observations and we need all the dims/etc to line up
    ans[[s]] <- list(zz = design@Z,
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

sparseToVec <- function(s, column = TRUE) {
  if (column) {
    as.matrix(s)[,1]
  } else {
    as.matrix(s)[1,]
  }
}
