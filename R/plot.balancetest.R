balanceTest_ggplot <- function(x,
                       xlab = "Standardized Differences",
                       absolute = FALSE,
                       strata.labels = NULL,
                       groups = NULL,
                       var.order = NULL,
                       var.grouping = NULL,
                       ...){
  
  if (!requireNamespace("ggplot2")) {
      stop("must have ggplot2 installed to use ggplot = TRUE plotting option")
  }
  # Tidyverse doesn't like rownames
  x <- tibble::as_tibble(x, rownames = "rowname", .name_repair = "unique")

  # Restructure data long
  tmp <- colnames(x)
  tmp <- tmp[tmp != "rowname"]
  x <- tidyr::pivot_longer(x, tmp, "strata")

  if (isTRUE(absolute)) {
    x <- dplyr::mutate(x, values = "abs(value)")
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

  if (!is.null(strata.labels)) {
    x$strata <- factor(x$strata, levels = strata.labels)
  }

  legend.title <- "Stratification" # just giving us some text for now
  plot <- ggplot2::ggplot(x,
                          ggplot2::aes_string(y = "rowname",
                                       x = "value",
                                       color = "strata",
                                       shape = "strata")) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::labs(color = legend.title,
                  shape = legend.title,
                  x = xlab,
                  y = ggplot2::element_blank())

    if (!is.null(x$group) && nlevels(x$group) > 0 ) {
        plot <- plot + ggplot2::facet_grid(group ~ ., scales = "free_y")
    }
  
  
  # Only draw connectors if there are multiple strata (avoids warning
  # from geom_path)
  if (length(unique(x$strata)) > 1) {
    plot <- plot + 
      ggplot2::geom_line(ggplot2::aes_string(group = "rowname"), color = "black")
  }
    
  return(plot)
}
