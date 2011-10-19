.wilcoxBackEnd <- function(
  adjusted.data,
  treatment, 
  blocks, 
  samples, p.value, summaries, # these will be ignored
  ...) {
  
  # TODO check blocks so see if they are paired
  # if (length(levels(blocks)) > 1) {
  #   strata <- blocks  
  # } else {
  #   strata <- NULL  
  # }
  
  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.data, 2, function(y) {
    res <- wilcox.test(y[!!treatment], y[!treatment])
    return(res[c("statistic", "p.value")])
  })

  tmp <- as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp) <- c("statistic", "p.value")

  return(new("RandomizationDistribution", 
      tmp, # RD inherits from data.frame
      test.statistic = wilcox.test,
      treatment = as.numeric(treatment),
      blocks = as.numeric(blocks)))
}

attr(wilcox.test, "randomizationEngine") <- .wilcoxBackEnd 
