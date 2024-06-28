#' Plot of balance across multiple strata
#'
#' The plot allows a quick visual comparison of the effect of different
#' stratification designs on the comparability of different
#' variables. This is not a replacement for the omnibus statistical test
#' reported as part of \code{\link{print.xbal}}. This plot does allow the
#' analyst an easy way to identify variables that might be the primary culprits
#' of overall imbalances and/or a way to assess whether certain important
#' covariates might be imbalanced even if the omnibus test reports that
#' the stratification overall produces balance.
#'
#' By default all variables and all strata are plotted. The scope
#' of the plot can be reduced by using the \code{\link{subset.xbal}} function to
#' make a smaller \code{xbal} object with only the desired variables or
#' strata.
#'
#' \code{\link{balanceTest}} can produce several different summary statistics for
#' each variable, any of which can serve as the data for this plot. By default,
#' the standardized differences between treated and control units makes a good
#' choice as all variables are on the same scale. Other statistics can be
#' selected using the \code{statistic} argument.
#'
#' The result of this function is a \code{\link{ggplot}} object. Most display of
#' the plot can be manipulated using additional commands appended to the plot
#' option. For example, the entire theme of the plot can be changed to black and
#' white using \code{plot(b) + theme_bw()}, where \code{b} is the result of a
#' call to \code{\link{balanceTest}}. The points on the plot are known as
#' "values", so colors or symbols used for each strata can be updated using the
#' \code{\link{scale_color_manual}} function. For example, \code{plot(b) +
#' scale_color_manaual(values = c('red', 'green', 'blue'))} for a balance test
#' of three stratification variables.
#'
#' @param x An object returned by \code{\link{xBalance}}
#' @param xlab The label for the x-axis of the plot
#' @param statistic The statistic to plot. The default choice of standardized
#' difference is a good choice as it will have roughly the same scale for all
#' plotted variables.
#' @param absolute Convert the results to the absolute value of the statistic.
#' @param strata.labels A named vector of the from \code{c(strata1 = "Strata Label 1", ...)}
#' that maps the stratification schemes to textual labels.
#' @param variable.labels A named vector of the from \code{c(var1 = "Var Label1", ...)}
#' that maps the variables to textual labels.
#' @param groups A vector of group names for each variable in
#' \code{x$results}. By default, factor level variables will be
#' grouped.
#' @param ... additional arguments to pass to \code{\link{balanceplot}}
#' @return A \code{ggplot2} object that can be further manipulated (e.g., to set the colors or text).
#' @seealso \code{\link{balanceTest}}, \code{\link{ggplot}}
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @import ggplot2
plot.balancetest <- function(x,
                             xlab = "Standardized Differences",
                             statistic = "std.diff",
                             absolute = FALSE,
                             strata.labels = NULL,
                             variable.labels = NULL,
                             groups = NULL,
                             ...) {
  ## this next line is for compatability with Josh's version of the function
  ## while I make it work with the usual plot signature
  var.order <- NULL; var.grouping <- NULL

  stopifnot(is.null(var.order) | is.null(var.grouping))

  tmp <- prepareXbalForPlot(x, statistic, absolute,
                            strata.labels, variable.labels) #square matrix, each row is row in plot, columns are horizontal locations
  #plot(tmp)
  autogroup <- attr(tmp, "term.labels")[attr(tmp, "groups")]
  autogroup[is.na(autogroup)] <- ""
  names(autogroup) <- rownames(tmp)
  
  if (is.null(groups)) {
    groups <- autogroup
  }
  x <- as.data.frame(tmp)
  balanceTest_ggplot(x, xlab = xlab, absolute = absolute, strata.labels = strata.labels, groups = groups)
}

balanceTest_ggplot <- function(x,
                       xlab = "Standardized Differences",
                       absolute = FALSE,
                       strata.labels = NULL,
                       groups = NULL,
                       var.order = NULL,
                       var.grouping = NULL,
                       ...) {

  # Tidyvers doesn't like rownames
  x <- tibble::as_tibble(x, rownames = "rowname", .name_repair = "unique")

  # Restructure data long
  x <- tidyr::pivot_longer(x, !"rowname", names_to = "strata")

  if (isTRUE(absolute)) {
    x <- dplyr::mutate(x, values = "abs(value)")
    if (xlab == "Standardized Differences") {
      xlab <- "Absolute Standardized Differences"
    }
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
      x$value <- as.numeric(x$value)
    }
  }

  x$group <- groups[as.character(x$rowname)]
  x$group[is.na(x$group)] <- ""
  x$group <- as.factor(x$group)

  if (!is.null(strata.labels)) {
    x$strata <- factor(x$strata, levels = strata.labels)
  }

  legend.title <- "Stratification" # just giving us some text for now
  plot <- ggplot2::ggplot(x,
                          ggplot2::aes(y = .data[["rowname"]],
                                       x = .data[["value"]],
                                       color = .data[["strata"]],
                                       shape = .data[["strata"]])) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::labs(color = legend.title,
                  shape = legend.title,
                  x = xlab,
                  y = ggplot2::element_blank())

    if (!is.null(x$group) && nlevels(x$group) > 0 ) {
      plot <- plot + ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y")
    }
  
  
  # Only draw connectors if there are multiple strata (avoids warning
  # from geom_path)
  if (length(unique(x$strata)) > 1) {
    plot <- plot +
      ggplot2::geom_line(ggplot2::aes(group = .data[["rowname"]]), color = "black")
  }
    
  return(plot)
}
