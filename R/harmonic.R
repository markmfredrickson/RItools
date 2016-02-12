##' Harmonic mean
##'
##' Calculate harmonic mean
##' @param data Data.
##' @return Harmonic mean.
harmonic <- function(data) {
  tapply(data$Tx.grp,data$stratum.code,function(x){2*sum( (x-mean(x))^2 )})
}
