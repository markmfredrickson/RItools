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

  predictions <- t(getApplyFunction(1:nrow(parameter.space), function(i) {
    do.call(invertModel, c(list(model, uniformity, z), parameter.space[i,]))
  }))

  std.parameters <- lapply(parameters, function(p) {
    minp <- min(p)
    maxp <- max(p)
    (p - minp) / (maxp - minp)
  })
  std.parameter.space <- do.call(expand.grid, std.parameters)

  # comparison order of dist() is
  # 1-2, 1-3, 1-4, ...,  2-3, 2-4, ... , 3-4, ...
  total.pairs <- dim(parameter.space)[1]

  left.index <- unlist(sapply(1:total.pairs, function(i) rep(i, total.pairs - i)))
  right.index <- unlist(sapply(2:total.pairs, function(i) i:total.pairs))

  pmat <- as.matrix(parameter.space)
  left <- pmat[left.index,, drop = F] ; colnames(left) <- paste("left", colnames(left), sep = ".")
  right <- pmat[right.index,, drop = F] ; colnames(right) <- paste("right", colnames(right), sep = ".")

  res <- cbind(left,
               right,
               parameter  = as.vector(dist(std.parameter.space)), 
               prediction = as.vector(dist(predictions)))

  class(res) <- c("parameterSensitivity", "matrix")
  attr(res, "parameters") <- parameters

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
  plot(hexbin(x[, c("parameter", "prediction")], ...))
}
