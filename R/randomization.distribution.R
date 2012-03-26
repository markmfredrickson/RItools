
###################################################
### Randomization distribution producing function
###################################################

# some utility classes
setClassUnion("OptionalList", c("list", "NULL"))
setClassUnion("OptionalDataFrame", c("data.frame", "NULL"))

# the RandomizationDistribution object
# for a single test statistic
# results of a call to randomizationDistribution
# each distribution is with respect to a given test statistic
setClass("RandomizationDistribution",
  representation(test.statistic = "function",
                 models.of.effect = "list", # first is always sharp null 
                 p.value = "function", # fn used to compute p-values
                 samples = "numeric", # the number of samples run, not necessarily requested
                 z  = "numeric", # 1/0 vector 
                 blocks = "numeric",
                 distribution = "matrix"), # if requested, the raw data
  contains = "data.frame")
# the matrix has a conventional form.
# The rows form the models tested, the columns are the samples
# The first column is the test statistic applied to the observed data.
# thus the number of samples is: dim(distrib)[2] - 1


randomizationDistributionEngine <- function(
  y,
  z,
  models, # a list of list(testStatistic, moe1, moe2, ...) all functions
  blocks = rep(1, length(z)),
  sampler = simpleRandomSampler(z = z, b = blocks),
  samples = 5000,
  p.value = general.two.sided.p.value,
  include.distribution = FALSE,
  summaries = list(),
  type = "exact",
  ...) {

  n <- length(z)
  if (is.null(blocks)) blocks <- rep(TRUE, n)
  blocks <- as.factor(blocks)
  stopifnot(n == length(z) && n == length(blocks))
 
  # note: it might be better to pass total pool, number of z, and a
  # vector of the size of the blocks. Rather than a vector indicating
  # z and a vector indicating block membership. keeping this signature
  # for now, but will think about it...
  randomizations <- sampler(samples)
  
  sharp.null <- function(y, z, b) { y }

  k <- length(models)

  apply.fn <- getApplyFunction()

  makedists <- function(i) {
    this.model <- models[[i]]
    test.statistic <- this.model[[1]]
    moes <- c(sharp.null, this.model[-1])
    
    # first, for each model, adjust the observed y with the observed
    # z
    adjusted.y <- sapply(moes, function(m) m(y, z, blocks, ...))

    # check if there is a backend function for this test.statistic
    if (type == "asymptotic" &&
        inherits(test.statistic, "AsymptoticTestStatistic")) {
      # basically, pass the adjusted y and everything else
      # to the backend and let it return a RandomizationDistribution
      # or subclass (probably a good idea)
      # backends may not honor samples, p.value, or summaries
      return(test.statistic@asymptotic(adjusted.y, z, blocks, 
                              samples, p.value, summaries, ...))  
    }
  
    # no backend, so use the standard approach
    # Possible todo: pull this out into its own "backend" function

    # now iterate over the randomizations, using the adjusted y
    this.distrib <- apply.fn(as.data.frame(randomizations$samples), function(z) {
      apply(adjusted.y, 2, function(d) { 
        test.statistic(d, z, blocks, ...)})})

    # see http://www.biostat.wustl.edu/archives/html/s-news/2000-01/msg00169.html
    # for a discussion of matrix + unlist vs. do.call + rbind
    this.distrib <- matrix(unlist(this.distrib), nrow = length(moes))

    adjusted.stats <- apply(adjusted.y, 2, function(d) { 
      test.statistic(d, z, blocks, ...)})

    pvs <- vector("numeric")

    for (i in 1:(length(moes))) { 
      pvs[i] <- p.value(adjusted.stats[i], this.distrib[i,])
    }
    
    this.result <- data.frame(statistic = adjusted.stats, p.value = pvs)

    if (length(summaries) > 0) {
      this.result <- cbind(this.result, lapply(summaries, function(f) { apply(this.distrib, 1, f)}))
    }
    
    tmp <- new("RandomizationDistribution", 
      this.result, # RD inherits from data.frame
      test.statistic = test.statistic,
      models.of.effect = moes,
      z = as.numeric(z),
      blocks = as.numeric(blocks),
      samples = dim(randomizations$samples)[2],
      p.value = p.value)

    if (include.distribution) {
      tmp@distribution <- this.distrib  
    }
    
    return(tmp)
  }

  distributions <- lapply(1:k, function(i){
      makedists(i)
  })
# temporary hack until I can get the intialize() method working
  names(distributions) <- names(models)

  return(distributions)
}

# front end as per Ben's recommendation. parameters are specified, not fns.
setClass("ParameterizedRandomizationDistribution",
  representation(
    call = "call",
    params = "data.frame"),
  contains = "RandomizationDistribution")

parameterizedRandomizationDistribution <- function(
  y,
  z,
  test.stat,
  moe = NULL, # single function with signature f(y, z, blocks, param1, param2, etc.)
  parameters = NULL, # list of name = values, name = values, ...
  ...) {
  
  # if either moe or parameters are present, both must be present
  if ((!is.null(moe) & is.null(parameters)) | (is.null(moe) & !is.null(parameters))) {
    stop("You must supply both parameters and a model effects if you supply either")
  }

  if (!is.null(moe) & !hasMethod(modelOfEffect, class(moe))) {
    stop("moe object must have a 'modelOfEffect' method (a model or a function)")  
  }

  if (!is.null(parameters) & !is.list(parameters)) {
    stop("Parameters must be a list")   
  }

  if (!is.null(parameters)) {
    parameter.space <- do.call(expand.grid, parameters)
    
    functions <- apply(parameter.space, 1, function(params) {
      force(params)
      function(y, z, blocks) {
        do.call(modelOfEffect, c(list(moe, y, z, blocks), params))
      }
    })

    rds <- randomizationDistributionEngine(y, z, models = list(c(test.stat,
      functions)), ...)
  } else {
    rds <- randomizationDistributionEngine(y, z, 
      models = list(c(test.stat)), ...)
  }

  rd <- as(rds[[1]], "ParameterizedRandomizationDistribution")

  if (!is.null(parameters)) {
    rd@params <- parameter.space
  }

  rd@call <- match.call()

  return(rd)
}

plot.ParameterizedRandomizationDistribution <- function(object, type = 'o', ...) {
  library(lattice)
  
  pnames <- colnames(object@params)
  np <- length(pnames)
  
  if (np == 0) {
    stop("Cannot plot sharp null only")
  }

  data <- cbind(object@params, object[-1,])

  if (np == 1) {
    fmla <- as.formula(paste("p.value ~ ", pnames[1]))
    return(xyplot(fmla, data = data, type = type, ...)) # drop the sharp.null
  }
  
  if (np == 2) {
    fmla <- as.formula(paste("p.value ~ ", pnames[1], "+", pnames[2]))
    return(levelplot(fmla, data = data, ...)) 
  }
  
  # this point, np > 2
  stop("Cannot plot parameters > 2 models (yet)")
}

setClass("ParameterizedRandomizationDistributionSummary",
  representation(
    randomizationDistribution = "ParameterizedRandomizationDistribution",
    point.estimate = "data.frame",
    showCall = "logical",
    sharp.null.p = "numeric"))

setMethod("summary", "ParameterizedRandomizationDistribution", function(object, showCall = T, ...) {
  return(new("ParameterizedRandomizationDistributionSummary", randomizationDistribution = object,
    point.estimate = point(object), showCall = showCall, sharp.null.p = object[1,2]))
})

setMethod("show", "ParameterizedRandomizationDistributionSummary", function(object) {
  if (object@showCall) {
    dc <- deparse(object@randomizationDistribution@call)
    cat("Call: ", dc[1], "\n")
    cat(paste("      ", dc[-1] , "\n", sep = ""))
    cat("\n")
  }

  tmp <- matrix(c(object@randomizationDistribution[1,1], object@sharp.null.p), 
    nrow = 1, byrow = T)

  rownames(tmp) <- c("Observed Test Statistic")
  colnames(tmp) <- c("Value", "Pr(>x)")
  printCoefmat(tmp, has.Pvalue = T, P.values = T) 
  cat("\n")

  # if there is more than one distribution (ie. more than the sharp null)
  # compute a point estimate
  if (dim(object@randomizationDistribution)[1] > 1) {
    cat("Hodges-Lehmann Point Estimate(s):\n")
    print(object@point.estimate)
  }

  invisible(object)
})

setGeneric("point",
  def = function(object) { standardGeneric("point") })
           
setMethod("point", "ParameterizedRandomizationDistribution", function(object) {

  combined <- cbind(object@params, object[-1,])
  # extract point estimate(s)
  maxp <- max(combined$p.value)
  point.estimate <- combined[combined$p.value == maxp, ]
  rownames(point.estimate) <- NULL
  return(point.estimate)
})


############################## Helper Functions ##############################

multicoreLoaded <- function() {
  "multicore" %in% loadedNamespaces()  
}

snowLoaded <- function() {
  "snow" %in% loadedNamespaces()  
}


getApplyFunction <- function() {
  opt <- options("RItools-apply")[[1]]

  if (!is.null(opt)) {
    return(opt)  
  }

  if (multicoreLoaded()) {
    return(mclapply) # yay speed!
  }

  if (snowLoaded()) { ## for now assumes that you start the cluster with cl<-makeCluster(...)
    return(function(x,fun){parLapply(cl,x,fun)}) # yay speed!
  }

  return(lapply) # the safe default
}

