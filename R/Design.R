################################################################################
# Design Objects: communicate cluster, strata, treatment, and covariate information
################################################################################

setClass("Design",
         representation = list(
             Z                 = "logical",
             Strata            = "data.frame",
             Cluster           = "factor",
             OriginalVariables = "character",
             Covariates        = "matrix"))
             

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
makeDesign <- function(fmla, data) {
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
    tbl <- table(str.data[, c(treatmentCol, s)])
    isGood <- apply(tbl, 2, function(x) { sum(x != 0) == 2 })
    if (!(all(isGood))) {
      stop("In ", s, ", the following stratum levels did not include both treated and control units: ",
           paste(collapse = ", ", colnames(tbl)[!isGood]))
    }
  }

  ## OK! data looks good. Let's proceed to make the design object with covariate data
  
  warning("TODO: add NA imputation when generating the covariates in makeDesign")
  data.fmla <- update(ts, paste("~", paste0(collapse = " - ", c(".", "1", vnames[str.idx]))))
  data.data <- model.frame(data.fmla, data, na.action = na.pass) #

  # we want our own contrast function for the factors that expands each level to its own dummy
  
  tlbl <- names(data.data)
  names(tlbl) <- as.character(tlbl)
  clist <- lapply(data.data, function(x) {
    if (is.factor(x)) { 
      structure(diag(nlevels(x)), dimnames = list(levels(x), levels(x)))
    } else {
      NULL
    }
  })
  clist <- clist[!sapply(clist, is.null)]
  
  data.mm   <- model.matrix(data.fmla, data.data, contrasts.arg = clist)

  # if we have clusters, aggegate up the counts and data, otherwise just return it as is.
  if (length(clusterCol) > 0) { 
    Cluster <- str.data[, clusterCol]
  } else {
    Cluster <- rep(1, dim(str.data)[1])
  }

  Z <- str.data[, treatmentCol]
  tmp <- str.data[, strataCols, drop = FALSE]
  colnames(tmp) <- gsub(colnames(tmp), pattern = "strata\\((.*)\\)", replacement = "\\1")
  Strata <- as.data.frame(lapply(tmp, factor))
  
  Covariates <- data.mm

  # create a look up table linking the model.matrix variables with the original variables
  
  originals <- attr(terms(data.fmla, data = data.data), "term.labels")[attr(data.mm, "assign")]

  return(new("Design",
             Z = as.logical(as.numeric(Z) - 1),
             Strata = Strata,
             Cluster = factor(Cluster),
             Covariates = Covariates,
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
                               design@Covariates)),
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
      hwts <- harmonic(data.frame(Tx.grp=zz,
                                  stratum.code=factor(design@Strata),
                                  data))
    }
    hwts <- hwts/sum(hwts, na.rm=TRUE)

    wtratio <- unsplit(sweights/hwts, design@Strata[[nn]], drop=TRUE)
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
weightedDesignToDescriptives <- function(design, covariate.scaling = TRUE) {
  stopifnot(inherits(design, "WeightedDesign")) # defensive programming

  res    <- list()
  vars   <- colnames(design@Covariates)
  strata <- colnames(design@Strata)

  s.p <- if (covariate.scaling) {
    xBalance.makepooledsd(design@Z, design@Covariates, length(vars))
  } else 1

  ans <- array(NA,
               dim = c(length(vars), 5, length(strata)),
               dimnames = list(
                   "vars" = vars,
                   "stat" = c("Control", "Treatment", "std.diff", "adj.diff", "adj.diff.null.sd"),
                   "strata" = strata))

  for (s in strata) {
    
    ss  <- design@Strata[, s]
    swt <- design@Weights[[s]]

    retain  <- !is.na(ss)
    ss      <- ss[retain]
    zz      <- design@Z[retain]
    covs    <- design@Covariates[retain, , drop = FALSE]
    wtratio <- swt$wtratio[retain]
    
    ### Calculate post.difference
    ZtH <- unsplit(tapply(zz, ss, mean), ss) ##proportion treated (zz==1) in strata s.
    ssn <- drop(crossprod(zz - ZtH, covs * wtratio)) ##weighted sum of mm in treated (using centered treatment)
    wtsum <- sum(unlist(tapply(zz, ss, function(x){ var(x) * (length(x) - 1) }))) ## h=(m_t*m_c)/m

    post.diff <- ssn/(wtsum) ##diff of means

    ans[, 'adj.diff', s] <- post.diff
    ans[, 'std.diff', s] <- post.diff/s.p

    ### Calculate post.Tx.eq.0, post.Tx.eq.1 --- now called "the treatment var"=0 and =1
    postwt0 <- unsplit(swt$sweights / tapply(zz == 0, ss, sum),
                       ss[zz == 0], drop=TRUE)

    ans[, "Control", s]   <- apply(covs[zz == 0,, drop = FALSE] * postwt0, 2, sum)
    ans[, "Treatment", s] <- ans[, "Control", s] + post.diff

    msmn <- xBalance.make.stratum.mean.matrix(ss, covs)

    tmat <- (covs - msmn)
    ##dv is sample variance of treatment by stratum
    dv <- unsplit(tapply(zz, ss, var), ss)
    ssvar <- colSums(dv * wtratio^2 * tmat * tmat) ## for 1 column in  mm, sum(tmat*tmat)/(nrow(tmat)-1)==var(mm) and sum(dv*(mm-mean(mm))^2)=ssvar or wtsum*var(mm)

    ans[, 'adj.diff.null.sd', s] <- sqrt(ssvar*(1/wtsum)^2)
  }
  
  return(ans)
}

# Aggregated Design totals up all the covariates and adds a new covariate "N"
#
setClass("AggregatedDesign",
         representation = list("N" = "numeric"),
         contains = "Design")

aggregateDesign <- function(design) {
  n.clusters <- nlevels(design@Cluster)
  Covariates <- matrix(NA, nrow = n.clusters, ncol = ncol(design@Covariates) + 1) # the extra col will be for cluster counts

  dupes <- duplicated(design@Cluster)
  Z <- design@Z[!dupes]
  Strata <- design@Strata[!dupes,, drop = FALSE]
  Cluster <- design@Cluster[!dupes]

  for (i in 1:length(Cluster))  {
    cname <- Cluster[i]
    subcovs <- design@Covariates[design@Cluster == cname, , drop = FALSE]
    Covariates[i, ] <- c(dim(subcovs)[1], colSums(subcovs))
  }

  colnames(Covariates) <- c("cluster.size", colnames(design@Covariates))
  
  new("AggregatedDesign",
      N = as.integer(Covariates[,1]),
      Z = Z,
      Strata = Strata,
      Cluster = Cluster,
      Covariates = Covariates,
      OriginalVariables = c("Cluster Size", design@OriginalVariables))
                       
}

alignDesignByStrata <- function(design, post.align.transform = NULL) {

  stopifnot(inherits(design, "WeightedDesign")) # defensive programming

  vars   <- colnames(design@Covariates)
  strata <- colnames(design@Strata)

  ans <- list()

  for (s in strata) {
    
    ss  <- design@Strata[, s]
    swt <- design@Weights[[s]]

    retain  <- !is.na(ss)
    ss      <- ss[retain]
    zz      <- design@Z[retain]
    covs    <- design@Covariates[retain, , drop = FALSE]
    wtratio <- swt$wtratio[retain]
    
    ### Calculate post.difference
    ZtH <- unsplit(tapply(zz, ss, mean), ss) ##proportion treated (zz==1) in strata s.
    ssn <- drop(crossprod(zz - ZtH, covs * wtratio)) ##weighted sum of mm in treated (using centered treatment)
    wtsum <- sum(unlist(tapply(zz, ss, function(x){ var(x) * (length(x) - 1) }))) ## h=(m_t*m_c)/m
    msmn <- xBalance.make.stratum.mean.matrix(ss, covs)

    tmat <- (covs - msmn)

    # dv is sample variance of treatment by stratum
    dv <- unsplit(tapply(zz,ss,var),ss)

    if (!is.null(post.align.trans)) {
      # Transform the columns of tmat using the function in post.align.trans
      tmat.new <- apply(tmat, 2, post.align.trans)

      # Ensure that post.align.trans wasn't something that changes the size of tmat (e.g. mean).
      # It would crash later anyway, but this is more informative
      if (is.null(dim(tmat.new)) || !all(dim(tmat) == dim(tmat.new))) {
        stop("Invalid post.alignment.transform given")
      }
      ## recenter on stratum means
      tmat <- tmat.new
      msmn <- xBalance.make.stratum.mean.matrix(ss, tmat)
      tmat <- tmat - msmn
      tmat <- tmat *swt$wtratio

      # Recalculate on transformed data the difference of treated sums and their null expectations
      # (NB: since tmat has just been recentered,
      # crossprod(zz,tmat) is the same as crossprod(zz-ZtH,tmat))
      ssn <- drop(crossprod(zz, tmat))
      ssvar <- apply(dv*tmat*tmat, 2, sum) 
    } else {
      tmat <- tmat *swt$wtratio
    }

    ans[[s]] <- tmat
  }
  return(ans)
}
