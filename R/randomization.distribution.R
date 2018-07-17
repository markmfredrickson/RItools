###################################################
### Randomization distribution producing function
###################################################

# some utility classes
setClassUnion("OptionalList", c("list", "NULL"))
setClassUnion("OptionalDataFrame", c("data.frame", "NULL"))

##' SharpNullTest class
##'
##' @slot observed.statistic observed.statistic
##' @slot samples sample
##' @slot call call
setClass("SharpNullTest",
  representation(observed.statistic = "numeric",
                 samples = "numeric",
                 call = "call"),
  contains = "array")

##' the RandomizationDistribution object
##'
##' for a single test statistic, results of a call to
##' randomizationDistribution, each distribution is with respect to a
##' given test statistic
##' @slot sharp.null SharpNullTest object
setClass("ParameterizedTest",
  representation(sharp.null = "SharpNullTest"),
  contains = "SharpNullTest")

##' Test a sharp hypothesis
##'
##' Test a sharp null hypothesis based on a experimental design.
##' Computes a test statistic of choice and reports p-values using direct
##' permutation
##' @param y Outcome variable of interest
##' @param z Treatment indicator vector, must be binary or logical
##' @param test.stat Test statistic to be reported
##' @param moe Model of effect, single function with signature 
##' f(y, z, param1, param2, etc.)
##' @param parameters Parameters values in model of effect 
##' List of name = values, name = values, ...
##' @param sampler Function specifying sampling procedure, uses simple 
##' randomization by default
##' @param samples Number of samples to use in direct permutation, uses 5000 
##' by default
##' @param type Either \code{"exact"} (direct permutation) or 
##' \code{"asymptotic"} (distributional approximation). Currently, 
##' \code{type = "asymptotic"} is not supported
##' @param ... Additional arguments
##' @return Observed test statistic and associated p-value
##' @seealso Update with vignettes here
##' @examples 
##' data(nuclearplants)
##' RItest(nuclearplants$cost, nuclearplants$pr, test.stat = t.mean.difference)
##' @export
RItest <- function(
  y,
  z,
  test.stat,
  moe = NULL,
  parameters = NULL,
  sampler = simpleRandomSampler(z = z, b = rep(1, length(z))),
  samples = 5000,
  type = "exact",
  ...) {

  # if either moe or parameters are present, both must be present
  if ((!is.null(moe) & is.null(parameters)) | (is.null(moe) & !is.null(parameters))) {
    stop("You must supply both parameters and a model effects if you supply either")
  }

  if (!is.null(moe) & !inherits(moe, "function")) {
    stop("moe object must be a function")
  }

  if (!is.null(parameters) & !is.list(parameters)) {
    stop("Parameters must be a list")
  }

  if (!(type %in% c("exact", "asymptotic"))) {
    stop("'type' argument must be either 'asymptotic' or 'exact'")
  }


  n <- length(z)

  # if(type == "exact") {
    randomizations <- sampler(samples) # a list with $weight and $samples args
  # }

  apply.fn <- getLApplyFunction

  # getting the p-value
  doit <- function(m) {
    function(...) {
      # first, for each model, adjust the observed y with the observed
      # z
      adjusted.y <- m(y, z, ...)
      obs.t <- test.stat(adjusted.y, z)

      # check if there is a backend function for this test.statistic
      # if (type == "asymptotic") {
      #     if (inherits(test.stat, "AsymptoticTestStatistic")) {
      #       # basically, pass the adjusted y and everything else
      #       # to the backend and let it return a RandomizationDistribution
      #       # or subclass (probably a good idea)
      #       # backends may not honor samples, p.value, or summaries
      #       return(test.stat@asymptotic(adjusted.y, z))
      #     }

      #   stop("No asymptotic backend exists for the test statistic")
      # }

      # no backend, so use the standard approach
      # Possible todo: pull this out into its own "backend" function

      # now iterate over the randomizations, using the adjusted y
      this.distrib <- apply.fn(as.data.frame(randomizations$samples), function(z) {
          test.stat(adjusted.y, z)})

      return(upper.p.value(obs.t, this.distrib))
    }
  }

  sharp.null <- function(y, z) { y }

  thecl <- match.call()

  sharp.null.test <- new("SharpNullTest",
                        doit(sharp.null)(),
                        observed.statistic = test.stat(y, z),
                        samples = dim(randomizations$samples)[2],
                        call = thecl)

  if (!is.null(moe)) {

    return(new("ParameterizedTest",
      farray(doit(moe), parameters), # inherits from array
      observed.statistic = sharp.null.test@observed.statistic,
      sharp.null = sharp.null.test,
      samples = sharp.null.test@samples,
      call = thecl))

  } else {
    return(sharp.null.test)
  }
}

#' Plot the results of RItest.
#'
#' @param x The results of RItest.
#' @param type The type (points, lines, etc) for one dimensional
#'   models. Ignored for multi-dimensional models.
#' @param summary A function to summarize higher dimension models (p >
#'   2) into 2 dimensions.
#' @param ... Add'l arguments
#' @return A lattice object. Use print to plot.
#' @export
#' @import lattice
#' @method plot ParameterizedTest
plot.ParameterizedTest <- function(x, type = 'o', summary = max, ...) {

  pnames <- dimnames(x)
  np <- length(pnames)

  if (np == 0) {
    stop("Cannot plot sharp null only")
  }

  if (np == 1) {
    flat <- data.frame(as.numeric(x), names(x))
    colnames(flat) <- c("p.value", pnames[1])
    fmla <- as.formula(paste("p.value ~ ", pnames[1]))
    return(xyplot(fmla, data = flat, type = type, ...)) # drop the sharp.null
  }

  return(levelplot(as.matrix(apply(x, 1:2, summary), ...)))

}

##' Show for SharpNullTest class
##'
##' @param object SharpNullTest
##' @return print call
##' @export
setMethod("show", "SharpNullTest", function(object) {
  print(object, showCall = TRUE)
})

##' Print for SharpNullTest
##'
##' @param x SharpNullTest
##' @param showCall Show call?
##' @param ... Ignored
##' @return printed object
##' @export
##' @method print SharpNullTest
print.SharpNullTest <- function(x, showCall = TRUE, ...) {
  if (showCall) {
    dc <- deparse(x@call)
    cat("Call: ", dc[1], "\n")
    cat(paste("      ", dc[-1] , "\n", sep = ""))
    cat("\n")
  }

  tmp <- matrix(c(x@observed.statistic, as.numeric(x)),
    nrow = 1, byrow = T)

  rownames(tmp) <- c("Observed Test Statistic")
  colnames(tmp) <- c("Value", "Pr(>x)")
  printCoefmat(tmp, has.Pvalue = T, P.values = T)
  cat("\n")

  invisible(x)
}

##' Print ParameterizedTest
##'
##' @param x A ParameterizedTest.
##' @param showCall Should call be shown?
##' @param ... Ignored
##' @return printed object
##' @export
##' @method print ParameterizedTest
print.ParameterizedTest <- function(x, showCall = TRUE, ...) {
  print.SharpNullTest(x@sharp.null, showCall)

  # compute a point estimate
  cat("Hodges-Lehmann Point Estimate(s):\n")
  print(point(x))
}

##' Generic point method
##'
##' @param object Object
##' @return something
##' @export
setGeneric("point",
           def = function(object) { standardGeneric("point") })

##' point for ParameterizedTest
##'
##' @param object ParameterizedTest
##' @return output
##' @export
setMethod("point", "ParameterizedTest", function(object) {

  # extract point estimate(s)
  maxp <- max(object)
  point.estimate.ind <- arrayInd(which(object == maxp),
                                 .dim = dim(object))

  point.estimate <- matrix(0,
                      nrow = dim(point.estimate.ind)[1],
                      ncol = dim(point.estimate.ind)[2])

  pset <- dimnames(object)

  for (i in 1:length(pset)) {
    point.estimate[, i] <- pset[[i]][point.estimate.ind[, i]]
  }

  colnames(point.estimate) <- names(dimnames(object))
  return(point.estimate)
})


############################## Helper Functions ##############################

##' Helper function to check if parallel package is loaded
##'
##' @return Logical
parallelLoaded <- function() {
   ("parallel" %in% loadedNamespaces())
}

##' Helper function to check if snow package is loaded
##'
##' @return Logical
snowLoaded <- function() { ## if a cluster named cl has been started
  length(find("cl"))==1 ##"snow" %in% loadedNamespaces()
}

##' Get proper lapply function
##'
##' @param x x
##' @param fun fun
##' @param ... Add'l arguments.
##' @return function
##' @import parallel
getLApplyFunction <- function(x,fun,...) {
  ## A parallized lapply
  opt <- getOption("RItools-lapply", lapply)

  if (!is.null(opt)) {
    return(opt(x,fun,...))
  }

  if (snowLoaded()) { ## for now assumes that you start the cluster with cl<-makeCluster(...)
    return(parLapply(cl,x,fun,...)) # yay speed!
  }

  if (parallelLoaded()) {
    options("mc.cores" = detectCores())
    return(mclapply(x,fun,...))
  }


  return(lapply(x,fun,...)) # the safe default
}


##' Get proper apply function
##'
##' @param x x
##' @param fun fun
##' @param ... Add'l arguments.
##' @return function
getApplyFunction <- function(x,fun,...) {
  ## A parallelized sapply
  opt <- getOption("RItools-sapply", sapply)

  if (!is.null(opt)) {
    return(opt(x,fun,...))
  }

  if (snowLoaded()) { ## for now assumes that you start the cluster with cl<-makeCluster(...)
    return(simplify2array(parLapply(cl,x,fun,...)))
  }

  if (parallelLoaded()) {
    options("mc.cores" = detectCores())
    return(simplify2array(mclapply(x,fun,...)))
  }


  return(sapply(x,fun,...)) # the safe default
}
