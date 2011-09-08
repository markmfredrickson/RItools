
###################################################
### Randomization distribution producing function
###################################################


produceRandomizations <- function(observed.treatment, blocks, samples) {
  if(!is.factor(blocks)){blocks<-factor(blocks)}
  # block size and the number of treated per block
  block.size <- table(blocks)
  block.treated <- aggregate(observed.treatment, list(blocks), sum)$x
  # these offsets will be needed later to turn relative into absolute indices
  block.starts <- append(0, cumsum(block.size))
  block.starts <- block.starts[1:(length(block.starts) - 1)]

  # overall statistics
  total.treated <- sum(observed.treatment)
  total.randomizations <- sum(lchoose(block.size, block.treated))

  # randomizations is a matrix (often abbreviated omega) of 
  # possible randomziations, for now ignoring blocks
  # it is generated either by direct enumeration (if small enough)
  # or by drawing from the distribution of randomizations
  if (total.randomizations > log(samples)) {
    randomizations <- matrix(nrow = total.treated, ncol = samples)

    for (i in 1:samples) {
      raw.draws <- mapply(sample.int, block.size, block.treated, SIMPLIFY = F)
      randomizations[, i] <- unlist(mapply(function(b,o) { b + o }, raw.draws, block.starts))
    NULL }
  } else { # small enough to figured exactly 
    # we first get all the block level samples (these vary in size
    block.combinations <- mapply(combn, block.size, block.treated, SIMPLIFY = F)
        
    block.combinations <- mapply(function(b,o) { b + o }, block.combinations, block.starts, SIMPLIFY = FALSE)
    # then all the indexes we need to get unique permutations of block.combs
    expansions <- 
      expand.grid(
        mapply(seq, 
               rep(1, length(block.size)), 
               sapply(block.combinations, function(b) { dim(b)[2] }),
               SIMPLIFY = F))
    
    # loop thru, creating a unique set of combinations
    exp.count <- dim(expansions)[1]

    randomizations <- matrix(nrow = total.treated, ncol = exp.count)
    for(i in 1:exp.count) {
      # from each block, grab the combination indexed by the row in expansions
      # combine all such items into a row in randomizations
      randomizations[,i] <- unlist(
        mapply(function(block, i) { block[, i] }, 
               block.combinations, 
               expansions[i,]))
      NULL
    }
  }

  return(randomizations)
}

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
                 treatment = "numeric", # 1/0 vector 
                 blocks = "numeric"),
  contains = "matrix")
# the matrix has a conventional form.
# The rows form the models tested, the columns are the samples
# The first column is the test statistic applied to the observed data.
# thus the number of samples is: dim(distrib)[2] - 1


randomizationDistributionEngine <- function(
  data,
  treatment,
  models, # a list of list(testStatistic, moe1, moe2, ...) all functions
  blocks = NULL,
  samples = 5000,
  ...) {

  n <- length(treatment)
  if (is.null(blocks)) blocks <- rep(TRUE, n)
  blocks <- as.factor(blocks)
  stopifnot(n == length(treatment) && n == length(blocks))
  
  # note: it might be better to pass total pool, number of treatment, and a
  # vector of the size of the blocks. Rather than a vector indicating
  # treatment and a vector indicating block membership. keeping this signature
  # for now, but will think about it...
  randomizations <- as.data.frame(produceRandomizations(treatment, blocks, samples))
  
  expand.z <- function(i) { a <- numeric(n); a[i] <- 1; return(a) }

  sharp.null <- function(y, z, b) { data }

  k <- length(models)

  apply.fn <- getApplyFunction()

  makedists <- function(i) {
    this.model <- models[[i]]
    test.statistic <- this.model[[1]]
    moes <- c(sharp.null, this.model[-1])
    
    # first, for each model, adjust the observed data with the observed
    # treatment indicator
    adjusted.data <- sapply(moes, function(m) m(data, treatment, blocks, ...))

    # now iterate over the randomizations, using the adjusted data
    this.distrib <- apply.fn(randomizations, function(z) {
      z <- expand.z(z)
      apply(adjusted.data, 2, function(d) { 
        test.statistic(d, z, blocks, ...)})})

    # see http://www.biostat.wustl.edu/archives/html/s-news/2000-01/msg00169.html
    # for a discussion of matrix + unlist vs. do.call + rbind
    this.distrib <- matrix(unlist(this.distrib), nrow = length(moes))

    # the format of a distribution is an m x k + 1 matrix, where
    # m is the number of models tested, and k is the number of randomizations
    # the extra column is the first column: test statistics under adjustment
    # to be used by the p-value functions
    adjusted.stats <- apply(adjusted.data, 2, function(d) { 
      test.statistic(d, treatment, blocks, ...)})
    this.distrib <- cbind(statistics = adjusted.stats, this.distrib)
    
    return(new("RandomizationDistribution", this.distrib,
      test.statistic = test.statistic,
      models.of.effect = moes,
      treatment = as.numeric(treatment),
      blocks = as.numeric(blocks)
      ))
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
  data,
  treatment,
  test.stat,
  moe = NULL, # single function with signature f(data, z, blocks, param1, param2, etc.)
  parameters = NULL, # list of name = values, name = values, ...
  blocks = NULL,
  samples = 5000) {
  
  # if either moe or parameters are present, both must be present
  if ((!is.null(moe) & is.null(parameters)) | (is.null(moe) & !is.null(parameters))) {
    stop("You must supply both parameters and a model effects if you supply either")
  }

  if (!is.null(moe) & !is.function(moe)) {
    stop("moe must be a function")  
  }

  if (!is.null(parameters) & !is.list(parameters)) {
    stop("Parameters must be a list")   
  }

  if (!is.null(parameters)) {
    parameter.space <- do.call(expand.grid, parameters)
    
    functions <- apply(parameter.space, 1, function(params) {
      force(params)
      function(data, z, blocks) {
        do.call(moe, c(list(data, z, blocks), params))
      }
    })

    rds <- randomizationDistributionEngine(data, treatment, models = list(c(test.stat,
      functions)), blocks, samples)
  } else {
    rds <- randomizationDistributionEngine(data, treatment, 
      models = list(c(test.stat)), blocks, samples)
  }

  rd <- as(rds[[1]], "ParameterizedRandomizationDistribution")

  if (!is.null(parameters)) {
    rd@params <- parameter.space
  }

  rd@call <- match.call()

  return(rd)
}

setClass("ParameterizedRandomizationDistributionSummary",
  representation(
    randomizationDistribution = "ParameterizedRandomizationDistribution",
    point.estimate = "data.frame",
    showCall = "logical",
    sharp.null.p = "numeric"))

setMethod("summary", "ParameterizedRandomizationDistribution", function(object, 
  p.value.function = general.two.sided.p.value, showCall = T, ...) {

  # point estimate(s)
  pvs <- p.values(object, p.value.function)
  maxp <- max(pvs$p)
  point.estimate <- pvs[pvs$p == maxp, ]
  rownames(point.estimate) <- NULL
  
  # observed test statistic is in the PRD object but we compute p-value against
  # the sharp null here
  sharp.null.p <- p.value.function(object[1,1], object[1,-1])

  return(new("ParameterizedRandomizationDistributionSummary", randomizationDistribution = object,
    point.estimate = point.estimate, showCall = showCall, sharp.null.p = sharp.null.p))
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

point.estimate.fn <- function(object,p.value.function = general.two.sided.p.value,...) {
  # extract point estimate(s)
  pvs <- p.values(object, p.value.function)
  maxp <- max(pvs$p)
  point.estimate <- pvs[pvs$p == maxp, ]
  rownames(point.estimate) <- NULL
  return(point.estimate)
}

setGeneric("point",def=function(object,p.value.function = general.two.sided.p.value,...) { standardGeneric("point") })
           
setMethod("point", "ParameterizedRandomizationDistribution",   point.estimate.fn)

setMethod("show", "ParameterizedRandomizationDistribution", function(object) {
  show(summary(object))
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

######################## xBalance Based Tests #################################


parameterizedXBalanceTest <- function(
  data,
  treatment,
  moe = NULL, # single function with signature f(data, z, blocks, param1, param2, etc.)
  parameters = NULL, # list of name = values, name = values, ...
  blocks = NULL,
  samples = 5000) {
  
  # if either moe or parameters are present, both must be present
  if ((!is.null(moe) & is.null(parameters)) | (is.null(moe) & !is.null(parameters))) {
    stop("You must supply both parameters and a model effects if you supply either")
  }

  if (!is.null(moe) & !is.function(moe)) {
    stop("moe must be a function")  
  }

  if (!is.null(parameters) & !is.list(parameters)) {
    stop("Parameters must be a list")   
  }

  if (!is.null(parameters)) {
    parameter.space <- do.call(expand.grid, parameters)
    
    functions <- apply(parameter.space, 1, function(params) {
      force(params)
      function(data, z, blocks) {
        do.call(moe, c(list(data, z, blocks), params))
      }
    })
  
    adjusted.data <- lapply(functions, 
      function(f) { f(data, treatment, blocks) })

    df <- data.frame(Z = treatment, unadjusted = data, adjusted.data)

    xb <- xBalance(Z ~ ., data = df, report = c("z.scores", "p.values"))
    output <- cbind(rbind(NA, parameter.space), p.value = xb$results[,,1][,2])

  } else {
    xb <- xBalance(Z ~ Y, data = data.frame(Z = treatment, Y = data), report =
      c("z.scores", "p.values"))
    output <- data.frame(NA, p.value = xb$results[,,1][,2])
  }

  return(output)
}
    

  

