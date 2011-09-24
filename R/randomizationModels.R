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
# direct.effect <- function(Z, blocks, ... parameters ...)
# spillover.effect <- function(Z, blocks, ... parameters ...)
#
# From these components, we can create modelsOfEffect (adjust observed data
# into the uniformity trial given a randomization and given parameters) and
# observedData (turning a uniformity trial into what we would expect given a
# randomization and some parameters).
#
################################################################################

### Models ###

setClass("RandomizationModel",
  representation(direct.effect = "function", 
                 spillover.effect = "function"))

randomizationModel <- function(direct.effect, spillover.effect) {
  new("RandomizationModel", direct.effect = direct.effect, spillover.effect = spillover.effect)  
}

# the NRM extends the RM by prepending the number of treated neighbors (ZS)
# onto the arguments of the direct.effect and spillover.effect functions
setClass("NetworkRandomizationModel",
  representation(S = "matrix"),
  contains = "RandomizationModel")

networkRandomizationModel <- function(S, direct.effect, spillover.effect) {
  de <- function(Z, blocks, ...) { direct.effect(Z, blocks, ZS = as.vector(Z %*% S), ...) }
  se <- function(Z, blocks, ...) { spillover.effect(Z, blocks, ZS = as.vector(Z %*% S), ...)}
  new("NetworkRandomizationModel", direct.effect = de, spillover.effect = se, S = S)
}


### Functions ###
setGeneric("modelOfEffect", function(model, R, Z, blocks = NULL, ...)
  standardGeneric("modelOfEffect"))

model.effect.helper <- function(model, combine, R, Z, blocks, ...) {
  tmp <- combine(R, Z * do.call(model@direct.effect, list(Z, blocks, ...)))
  combine(tmp, (1 - Z) * do.call(model@spillover.effect, list(Z, blocks, ...))) 
}

setMethod("modelOfEffect", signature = c("RandomizationModel"),
function(model, R, Z, blocks = NULL, ...) {
  model.effect.helper(model, combine = `-`, R, Z, blocks, ...)
})

setGeneric("observedData", function(model, R, Z, blocks = NULL, ...)
  standardGeneric("observedData"))

setMethod("observedData", signature = c("RandomizationModel"),
function(model, R, Z, blocks = NULL, ...) {
  model.effect.helper(model, combine = `+`, R, Z, blocks, ...)
})

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

# Yt = Yc + Z * tau
constant.additive.model <- randomizationModel(returnTau, returnZero)

## The following models are parameterized with a network matrix S

# the linear spillover model from Bowers and Fredrickson, 2011
# Yt = Yc + Z * tau + (1 - Z) * min(tau, beta * Z^t %*% S)
linear.spillover.fn <- function(Z, blocks, ZS, tau, beta) {
  pmin(tau, beta * ZS) 
}

linear.spillover.model <- function(S) {
  networkRandomizationModel(S, returnTau, linear.spillover.fn) 
}

# the inverse growth spillover model from Bowers and Fredrickson, 2011
inverse.spillover.fn <- function(Z, blocks, ZS, tau1, tau2) {
  # Growth curve: when tau2<0 effect is decreasing in the number of treated
  # neighbors, 
  # when tau2>0, effect is increasing in the number of treated neighbors, 
  # tau2=0 means all same effect regardless of the neighborhood. 
  # Limited to  effect no greater than tau1
  tau1 * (1 - (1 / (1 + ZS ^(tau2))))
}

inverse.spillover.model <- function(S) {
  networkRandomizationModel(S, returnTau, inverse.spillover.fn)
}


### Analyzing Operating Characteristics ###


