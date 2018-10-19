##' Harmonic mean
##'
##' Calculate harmonic mean
##' @param data Data.
##' @return numeric vector of length \code{nlevels(data$stratum.code)}
##' @keywords internal
harmonic <- function(data) {
  tapply(data$Tx.grp,data$stratum.code,function(x){2*sum( (x-mean(x))^2 )})
}
##' Harmonic mean times mean of weights
##'
##' 
##' @param data 
##' @return numeric vector of length \code{nlevels(data$stratum.code)}
##' @keywords internal
harmonic_times_mean_weight <- function(data) {
    weightmeans <- tapply(data$unit.weights,data$stratum.code, mean)
    hmeans <- harmonic(data)
    stopifnot(all(names(hmeans)==names(weightmeans))) 

    weightmeans * hmeans
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

