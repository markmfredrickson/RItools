##' Harmonic mean
##'
##' Calculate harmonic mean
##' @param data Data.
##' @return Harmonic mean.
harmonic <- function(data) {
  tapply(data$Tx.grp,data$stratum.code,function(x){2*sum( (x-mean(x))^2 )})
}

effectOfTreatmentOnTreated <- function(data) {
  tapply(data$Tx.grp, data$stratum.code, function(x) {
    sum(x > 0) / length(x)
  })
}
