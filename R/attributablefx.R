################################################################################
# Attributable effects: Models of effect that apply to aggregate values
################################################################################

##' Attribute Effects helper functions
##'
##' @param y y
##' @param z z
##' @param b b
##' @param a a
##' @return value
##' @name attributable.effect.helpers
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

##' @rdname attributable.effect.helpers
max.attributable.effects <- function(y, z) {
  sum(y == 1 & z == 1)
}

##' @rdname attributable.effect.helpers
min.attributable.effects <- function(y, z) {
  0 - sum(y == 0 & z == 1)
}

##' Attributed Effects
##'
##' @param data data
##' @param treatment treatment
##' @param test.stat test.stat
##' @param samples samples
##' @param min.a min.a
##' @param max.a max.a
##' @param step step
##' @return RITest
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
