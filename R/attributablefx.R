################################################################################
# Attributable effects: Models of effect that apply to aggregate values
################################################################################

# the main model: assumes a binary outcome
additive.attributable.effect <- function(y, z, b, a) {
  if (a == 0) {
    return(y)
  }
  
  z <- as.logical(z)

  if (a > 0) {
    y[z & y == 1][1:a] <- 0
  } else {
    y[z & y == 0][1:(-a)] <- 1
  }

  return(y)                                                
}

# helper functions to compute min and max possible attributable effects
max.attributable.effects <- function(y, z) {
  sum(y == 1 & z == 1)
}

min.attributable.effects <- function(y, z) {
  0 - sum(y == 0 & z == 1)
}

attributableEffects <- function(
  data,
  treatment,
  test.stat,
  samples = 5000,
  min.a = min.attributable.effects(data, treatment),
  max.a = max.attributable.effects(data, treatment),
  step = 1) {

  if (!is.numeric(data) | !all(data %in% c(0,1))) {
    stop("Argument \'data\' must be numeric 0/1")
  }

  RItest(data, 
    treatment,
    test.stat,
    additive.attributable.effect,
    parameters = list(a = seq(min.a, max.a, step)),
    samples)
  
}
