################################################################################
# Smoothing the already smooth integration with C++ via Rcpp
################################################################################

#' @imports RcppExports.R
setClass("RcppFunction",
         representation = c("RcppFn" = "externalptr"),
         contains = "function")

# RcppExports.R has the definitions for the R callable functions.
# For the ones that might also get called as functions in C++ algorithms, 
# the objects are promoted the RcppFunction class with the extra slot.
# Since we have to get the pointer at run time, we need to wait for onLoad in order 
# actually populate the RcppFunction slot.

#' @imports testStatistics.R
meanDifference <- as(meanDifference, "RcppFunction")
wilcoxTestStatistic <- as(wilcoxTestStatistic, "RcppFunction")

#' @imports Rcpp.R
.onLoad <- function(...) {
  message("Loading")
  meanDifference@RcppFn <<- testStatisticPtr("meanDifference")
  wilcoxTestStatistic@RcppFn <<- testStatisticPtr("wilcoxTestStatistic")
}
