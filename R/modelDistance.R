#' Compares the predictions of a model over the parameter space.
#'
#' Models that are insensitive to changes in the parameters will have little distance between predictions.
#' 
#' @param model A \code{\link{UniformityModel}} object that will be used with the parameters to generate the model predictions.
#' @param parameters A list of the form \code{list(param1 = c(1,2,3), param2 = c(4,5,6), ...)} where the names are the parameter arguments to the model function.
#' @param uniformity Data that represent the uniformity trial (the state of the world when treatment is not applied to any unit).
#' @param z A vector of treatment indicators.
#' @return A matrix of two columns: the parameter space distance and the prediction distance for all pair-wise comparisons
#' @export
parameterSensitivity <- function(model, parameters, uniformity, z) {

  parameter.space <- do.call(expand.grid, parameters)

  predictions <- t(apply(parameter.space, 1, function(r) {
    do.call(invertModel, c(list(model, uniformity, z), r))
  }))

  res <- matrix(c(as.vector(dist(parameter.space)), as.vector(dist(predictions))), ncol = 2)
  colnames(res) <-  c("parameter", "prediction")

  class(res) <- c("parameterSensitivity", "matrix")

  return(res)
}

setOldClass(c("parameterSensitivity", "matrix"))

#' Plots the results of \code{\link{parameterSensitivity}} as a 2D histogram.
#'
#' @param x The matrix of distances from \code{\link{parameterSensitivity}}.
#' @param ... Arguments to be passed to the plot method for \code{\link{hexbin::hexbin}}
#' @import hexbin
#' @export
plot.parameterSensitivity <- function(x, ...) {
  plot(hexbin(x), ...)
}
