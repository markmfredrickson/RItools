################################################################################
# Models of Effect functions for Randomization Distributions
#
# A moe is a function of at least 3 arguments
# - ys: the data to be adjusted, these methods assume a vector
# - z: a numeric indicator of treament (0,1)
# - b: a blocking indicator (may be a single level)
#
# Additional arguments are considered population parameters of the model.
# See RItest for more details.
################################################################################



##' Helper function
##'
##' @param ys the data to be adjusted, these methods assume a vector
##' @param z a numeric indicator of treament (0,1)
##' @param b a blocking indicator (may be a single level)
##' @param beta beta
##' @return value
constant.multiplicative.model <- function(ys, z, b, beta) {
  ys / ((z * beta) + (1 - z))
}

# This function is not yet ready for prime time. The idea is that after adjusting
# the data, a regression is run on the _adjusted_ data. Then the residuals are returned
# for the test statistic function to consume.
# This should be a test statistic...
# lm.residualizer take a formula of the form y ~ ..., with all values
# included in the ... included in the data argument (a data.frame).
# it returns a function that takes new y values and applies the model to them.
# the inner function returns the residuals of the model for later test stats.
### lm.residual.moe <- function(model, fmla, data){
###   # it might be more efficient to build the model once, and then use
###   # update(). I'm proceeding with the simplier implementation for now
###   function(...) {
###     y <- residuals(lm(fmla, data))
###     model <- lm(fmla, cbind(y = adjusted.outcomes, data))
###     residuals(model)
###   }
### }


##' Min Max Model
##'
##' Given lower and lower bounds a and b, create a new model that
##' limits within those bounds
##' @param model Model
##' @param lower Lower bound
##' @param upper Upper bound
##' @return New model
min.max.model <- function(model, lower = -Inf, upper = Inf) {
  function(...) {
    pmin(upper, pmax(model(...), lower))
  }
}
