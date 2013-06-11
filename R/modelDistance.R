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

#' Plots the results of the parameter sensitivity on the original parameter space.
#' 
#' Works best for 2D models, may not work for other dimenions.
#'
#' @param x The results of the call to \code{\link{parameterSensitivity}}.
#' @param subset An expression that will evaluate to a logical value when
#' evaluated in the context of \code{x}. You can use \code{parameter} and
#' \code{prediction} to subset which parameters to show.
#' @export
parameterSensitivityBreakoutPlot <- function(x, subset) {
  params <- attr(x, "parameters")
 
  ps <- length(params)

  if (ps > 2) {
    stop("Cannot plot models with more than 3 parameters")
  }

  if (ps < 1) {
    stop("Cannot plot models with no parameters")
  }
  
  if (missing(subset)) 
    r <- TRUE
  else {
    e <- substitute(subset)
    r <- eval(e, as.data.frame(x), parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must be logical")
    r <- r & !is.na(r)
  }

  # select only the relevant portion of x and only the parameter columns
  x.subset <- x[r, 1:(2 * ps)] 
  
  left <- x.subset[, 1:ps]
  right <- x.subset[, (ps + 1):(2 * ps)]

  xlims <- range(params[[1]])
  ylims <- range(params[[2]])

  plot(NULL, xlim = xlims, ylim = ylims, xlab = names(params)[1], ylab = names(params)[2])
  segments(left[,1] , left[,2],
           right[,1], right[,2])
}
