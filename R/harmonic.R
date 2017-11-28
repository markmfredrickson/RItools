##' Harmonic mean
##'
##' Calculate harmonic mean
##' @param data Data.
##' @return Harmonic mean.
##' @keywords internal
harmonic <- function(data) {
  tapply(data$Tx.grp,data$stratum.code,function(x){2*sum( (x-mean(x))^2 )})
}

##' Number of treatment clusters by stratum
##' 
##' Calculate the number of treatment clusters by stratum --
##' these being proportional to "effect of treatment on treated" weights
##' when assignment probabilities are uniform within each stratum.
##' 
##' NB: currently, i.e. as of this inline note's commit,
##' this function is used only in testing.
##'
##' @param data Data.
##' @return Cluster count
##' @keywords internal
##'
effectOfTreatmentOnTreated <- function(data) {
    txcts <- table(as.logical(data$Tx.grp), data$stratum.code)
    txcts["TRUE",]
}

