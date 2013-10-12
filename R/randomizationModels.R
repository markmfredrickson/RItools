################################################################################
# Models for use in randomization distributions
# 
# All of our models take the same basic form:
#
# Observed Data - Z * f(treated neighbors, ... params ... ) - (1 - Z) *
# g(treated neighbors, ... params)
#
# Model objects can encapsulate this basic using just
# direct.effect and spillover.effect functions.
# 
# direct.effect <- function(Z, ... parameters ...)
# spillover.effect <- function(Z, ... parameters ...)
#
# From these components, we can create modelsOfEffect (adjust observed data
# into the uniformity trial given a randomization and given parameters) and
# observedData (turning a uniformity trial into what we would expect given a
# randomization and some parameters).
#
################################################################################

### Models ###


## The most generic type of model is a "UniformityModel" which is a function
## that takes y and z (and possibly parameters) and returns the uniformity trial
## the @inverse slot contains a function that given the uniformity trial data, z, and params
## returns the observed data implied by the model.
setClass("UniformityModel", 
  contains = "function",
  representation(inverse = "function"))

UniformityModel <- function(f, inv) { new("UniformityModel", f, inverse = inv)}

## apply the inverse of the model to the uniformity trial data to see what
## we would observe given z and parameters
invertModel <- function(m, y_0, z, ...) {
  m@inverse(y_0, z, ...)  
}

setGeneric("givenParams", function(model, ...)
  standardGeneric("givenParams"))

setMethod("givenParams", "function",
function(model, ...) {
  dots <- match.call(expand.dots = FALSE)[["..."]]
  function(y, z, b, ...) {
    do.call(model, c(list(y, z, b, ...), dots))
  }
})

setMethod("givenParams", "UniformityModel",
function(model, ...) {
  dots <- match.call(expand.dots = FALSE)[["..."]]
  UniformityModel(
    function(y, z, ...) {
      do.call(model, c(list(y, z, ...), dots))
    },
    function(y0, z, ...) {
      do.call(invertModel, c(list(model, y0, z), dots))  
    })
})

## levelEffectModel: for each level in Z, applies a specific function to create uniformity
## each function is itself a model.
## the goal here is to generalize out the y - z_1 * tau1 - z_2 * tau2 - ... idea for arbitrary fns
## setClass("LevelEffectModel",
##   contains = "UniformityModel",
##   representation(levelFunctions = "list"))
## 
## LevelEffectModel <- function(effects) {
##   if(!inherits(effects, "list") && !(all(sapply(effects, function(x) { inherits(x, "function")})))) {
##     stop("All effects must be functions")
##   }
##   
##   if (is.null(names(effects))) {
##     names(effects) <- 0:(length(effects) - 1)  
##   }
## 
##   nms <- names(effects)
## 
##   if (all(sapply(effects, function(x) { inherits(x) }))) {
##     return( UniformityModel(
##       function(y, z, ...) {
##         for (i in 1:length(effects)) {
##           lvl <- nms[i]
##           y <- y - effects[[i]](y, as.numeric(z == lvl), ...)  
##         }
##         y 
##       },
##       function(y_0, z, ...) {
##         for (i in 1:length(effects)) {
##           lvl <- nms[i]
##           y_0 <- y_0 + (z == lvl) * effects[[i]](y, ...)  
##         }
##         y_0
##       }))
##   }
##   
## }

# the NRM extends the RM by prepending the number of treated neighbors (ZS)
# onto the arguments of the direct.effect and spillover.effect functions
# setClass("NetworkRandomizationModel",
#   representation(S = "matrix"),
#   contains = "RandomizationModel")
# 
# networkRandomizationModel <- function(S, direct.effect, spillover.effect) {
#   de <- function(Z, ...) { direct.effect(Z, ZS = as.vector(Z %*% S), ...) }
#   se <- function(Z, ...) { spillover.effect(Z, ZS = as.vector(Z %*% S), ...)}
#   new("NetworkRandomizationModel", direct.effect = de, spillover.effect = se, S = S)
# }


model.effect.helper <- function(model, combine, R, Z, ...) {
  tmp <- combine(R, Z * do.call(model@direct.effect, list(Z, ...)))
  combine(tmp, (1 - Z) * do.call(model@spillover.effect, list(Z, ...))) 
}


# helper function for all the models that have a "tau" direct effect
returnTau <- function(...) {
  arguments <- match.call(expand.dots = TRUE)
  return(arguments[["tau"]])
}

returnZero <- function(...) {
  return(0)  
}

# I need to figure out how to handle Z, which is required for the direct.effect
# and spillover.effect models.
# 
# plot.RandomizationModel <- function(object, params, ...) {
#   library(lattice)
#   
#   pgrid <- do.call(expand.grid, params)
# 
#   primary <- pgrid[, 1]
#   primary.name <- colnames(pgrid)[1]
#   others <- pgrid[, -1]
#   other.names <- colnames(pgrid)[-1] 
# 
#   xlim <- range(primary)
#   
#   direct.lines <- apply(others, 1, function(row) {      
#     f <- function(x) { do.call(object@direct.effects, )}
#   })
#   
#   
# }

### Common Models ###

# World's simplest model, the sharp null of no effects
sharp.null.model <- UniformityModel(function(y, z) { y }, function(y, z) { y })

# Yt = Yc + Z * tau
constant.additive.model <- UniformityModel(function(y, z, tau) { y - z * tau }, 
                                           function(y_0, z, tau) { y_0 + z * tau})

## The following models are parameterized with a network matrix S
## network models are commented out for now
# the linear spillover model from Bowers and Fredrickson, 2011
# Yt = Yc + Z * tau + (1 - Z) * min(tau, beta * Z^t %*% S)
# linear.spillover.fn <- function(Z, ZS, tau, beta) {
#   pmin(tau, beta * ZS) 
# }
# 
# linear.spillover.model <- function(S) {
#   networkRandomizationModel(S, returnTau, linear.spillover.fn) 
# }
# 
# # the inverse growth spillover model from Bowers and Fredrickson, 2011
# inverse.spillover.fn <- function(Z, ZS, tau1, tau2) {
#   # Growth curve: when tau2<0 effect is decreasing in the number of treated
#   # neighbors, 
#   # when tau2>0, effect is increasing in the number of treated neighbors, 
#   # tau2=0 means all same effect regardless of the neighborhood. 
#   # Limited to  effect no greater than tau1
#   tau1 * (1 - (1 / (1 + ZS ^(tau2))))
# }
# 
# inverse.spillover.model <- function(S) {
#   networkRandomizationModel(S, returnTau, inverse.spillover.fn)
# }


### Analyzing Operating Characteristics ###
# powerAnalysis: 
#  model: a RandomizationModel object, e.g. constant.additive.model
#  Z: an example randomization that will be used to to generate possible
#    possible randomizations 
#  true.params: list(param.name = value, ...) representing the "true"
#    parameters of this simulation  
#  test.params: list(param.name = c(1,2,3...), ...) representing "false"
#    values to test for rejection
#  data: simulated data from a uniformity trial (i.e. all get control)
#  test.samples: number of samples to draw from the randomization distribution
#    when testing a hypothesis
#  power.samples: number of datasets to generate simulated Z values
#  ... arguments to be passed to pRD (and probably then to the engine)

## Commented out until there is a replacement for randomizationDistributionEngine
# compareModels <- function(models, # a list of (one param) models to try 
#                           repetitions, # how many times to create data
#                           test.statistic,                          
#                           uniformity, 
#                           sampler,
#                           summary.column = "p.value",
#                           ...) {
#   
#   Zs <- sampler(repetitions)
# 
#   if (missing(uniformity)) {
#     uniformity <- rnorm(length(Zs$samples[,1]))  
#   }
# 
#   results <- lapply(as.data.frame(Zs$samples), function(z) {
#     data <- lapply(models, function(m) { invertModel(m, uniformity, z)})
# 
#     test <- lapply(data, function(d) {
#       tmp <- randomizationDistributionEngine(d, z, list(c(test.statistic, models)),
#         sampler = sampler, ...) # all other args passed along, e.g. samples  
# 
#       tmp[[1]][-1, summary.column, drop = F] # drop the sharp null row
#     })
#     
#     test <- matrix(unlist(test), nrow = length(models), dimnames = list(names(models), names(models)))
#     return(test)
#   })
# 
#   names(results) <- 1:repetitions
# 
#   require(abind)
# 
#   return(abind(results, along = 3))
# }
# 
# parameterizedCompareModels <- function(models, ...) {
#   unparameterized <- lapply(models, function(l) {
#     m <- l[[1]]
#     params <- l[-1]
# 
#     if (length(params) > 0) {
#       parameter.space <- do.call(expand.grid, params)
#     return(apply(parameter.space, 1, function(ps) {
#         force(ps)
#         do.call(givenParams, c(list(m), ps))
#         }))
#     } else {
#       return(m)  
#     }
#   }) 
#   
#   compareModels(unlist(unparameterized), ...)
# }

# given the results of analyzeModel, what can we say?

testSize <- function(analysis, alpha = NULL, tol = .Machine$double.eps ^ 0.5) {
  # unpack the analysis object into some useful components
  true.params <- analysis$true.params
  sims <- analysis$simulation # a list of pRD objects
  
  smash <- sapply(sims, function(s) { s$p.value })

  # if alpha is not passed, use all unique p-values found in the
  # simulation
  if (is.null(alpha)) {
    alpha <- sort(unique(as.vector(smash)))
  }

  # what were the p-values listed for the true paramers?
  # the additional false is for the null model which is the first row of every pRD obj
  pspace <- sims[[1]]@params
  psnames <- names(pspace)
  idx <- c(FALSE, apply(pspace, 1,function(x){ identical(as.numeric(x), as.numeric(true.params[psnames]))}))

  p.truths <- smash[idx,]
   
  res <- t(sapply(alpha, function(a) { 
    c(nominal.alpha = a, realized.alpha = mean((p.truths - a) <= tol))
  }))

  return(res)

}
