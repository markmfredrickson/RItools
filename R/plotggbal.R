#' ggplot2 based balance plot
#'
#' @param x Balance object
#' @param xlab Same as plot.xbal
#' @param statistic Same as plot.xbal
#' @param colors Vector of colors, same length as number of strata
#' @param symbols  Vector of (numeric) symbols, same length as number of strata
#' @param absolute Same as plot.xbal
#' @param strata.labels Same as plot.xbal
#' @param variable.labels Same as plot.xbal
#' @param legend.position Where the legend should go ("bottom", "right")
#' @param legend.title Title of legend (blank by default)
#' @param var.order Vector of variable names (variable.labels if non-NULL) defining ordering. Otherwise alphabetical. Only one of `var.order` and `var.groupings` can be non-NULL.
#' @param var.groupings List of named vectors of variables names (variable.labels if non-NULL) which should be grouped. The name of the list entry is the grouping label. Only one of `var.order` and `var.groupings` can be non-NULL.
#' @param ...
#'
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @import ggplot2
plotggxbal <- function(x,
                       xlab = "Standardized Differences",
                       statistic = "std.diff",
                       colors = NULL,
                       symbols = NULL,
                       absolute = FALSE,
                       strata.labels = NULL,
                       variable.labels = NULL,
                       # groups = NULL, NYI
                       legend.position = "bottom",
                       legend.title = "",
                       var.order = NULL,
                       var.grouping = NULL,
                       ...) {

  stopifnot(is.null(var.order) | is.null(var.grouping))

  x <- as.data.frame(RItools:::prepareXbalForPlot(x, statistic, absolute,
                                        strata.labels, variable.labels))

  # Tidyverse doesn't like rownames
  x <- tibble::rownames_to_column(x)

  # Restructure data long
  x <- tidyr::gather(x, strata, values, -rowname)

  if (isTRUE(absolute)) {
    x <- dplyr::mutate(x, values = abs(values))
  }

  if (!is.null(var.order)) {
    stopifnot(length(var.order) == length(unique(x$rowname)))
    x$rowname <- factor(x$rowname,
                        levels = rev(var.order))
  } else {
    x$rowname <- factor(x$rowname,
                        levels = sort(unique(x$rowname),
                                      decreasing = TRUE))
  }

  if (!is.null(var.grouping)) {
    for (grouplabel in names(var.grouping)) {
      curorder <- levels(x$rowname)
      neworder <- curorder[!(curorder %in% var.grouping[[grouplabel]])]
      neworder <- c(rev(var.grouping[[grouplabel]]), grouplabel, neworder)
      x$rowname <- factor(x$rowname, levels = neworder)
      for (stratas in unique(x$strata)) {
        x <- rbind(x, c(grouplabel, stratas, NA))
      }
      x$values <- as.numeric(x$values)
    }
  }





  plot <- ggplot2::ggplot(x,
                          ggplot2::aes(y = rowname,
                                       x = values,
                                       color = strata,
                                       shape = strata)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_line(ggplot2::aes(group = rowname), color = "black") +
    ggplot2::geom_point() +
    ggplot2::theme(legend.position = legend.position) +
    ggplot2::labs(color = legend.title,
                  shape = legend.title,
                  x = xlab,
                  y = ggplot2::element_blank())

  if (!is.null(colors)) {
    stopifnot(length(colors) == length(unique(x$strata)))
    plot <- plot + ggplot2::scale_color_manual(values = colors)
  }
  # Should probably also take colors/symbols arguments of length 1 and expand them.
  if (!is.null(symbols)) {
    stopifnot(length(symbols) == length(unique(x$strata)))
    plot <- plot + ggplot2::scale_color_manual(values = symbols)
  }

  return(plot)
}
