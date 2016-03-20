##' Harmonic mean
##'
##' Calculate harmonic mean
##' @param data Data.
##' @return Harmonic mean.
harmonic <- function(data) {
  tapply(data$Tx.grp,data$stratum.code,function(x){2*sum( (x-mean(x))^2 )})
}

effectOfTreatmentOnTreated <- function(data) {
    txcts <- table(as.logical(data$Tx.grp), data$stratum.code)
    txcts["TRUE",]/(txcts["TRUE",] + txcts["FALSE",])
  }


