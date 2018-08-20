################################################################################
##
## Classes and functions related to things that can be assigned to treatment.
##
################################################################################


setClassUnion("OptionalDataFrame", c("NULL", "data.frame"))

## A virtual class that represents things on which treatment assignments can be made.
setClass("Assignable", contains = "OptionalDataFrame")


#' The total number of units in the study. May be larger than the number of
#' clusters that can be assigned to treatment.
#'
#' @param x The assignable object.
#' @return The number of individual units in the study. For subgroups, can be a vector of counts
setGeneric("number_of_units", def = function(x) { NA })

#' The total number of assignable objects that can be assigned to treatment, may
#' be smaller than the number of units.
#'
#' @param x The assignable object.
#' @return The number of objects that can be assigned to treatment.
setGeneric("number_of_assignments", def = function(x) { NA })

## Units represent the lowest possible grouping in a design. They can have data.
## For a set of units, we might only care about a subgroup of them at any given time.
## For now, we carry that information around with us.
setClass("Units", contains = "Assignable",
         slots = c(subgroup = "factor"))

## TODO: add a way to handle just an integer
units <- function(covariates, subgroup = NULL) {
    n <- dim(covariates)[1]
    if (is.null(subgroup)) {
        subgroup <- rep(1, n)
    }

    if (n != length(subgroup)) {
        stop("If a subgroup is included it must refer to all units.")
    }

    new("Units", covariates, subgroup = subgroup)
}

setMethod("number_of_units", signature = c("Units"), function(x) { sum(!is.na(x@subgroup)) })
setMethod("number_of_assignments", signature = c("Units"), function(x) { dim(x)[1] })

setClass("Clusters",
         contains = "Assignable",
         slots = list(map = "integer", ## TODO: should this be a factor instead?
                      covariates = "OptionalDataFrame"))

nclusters <- function(x) { length(unique(x)) }

#' Clustered treatment assignment mechanisms.
#'
#' Units are arranged into clusters, which form the basis of randomization.
#' Units are mapped to clusters when computing quantities. The clusters require
#' their own design object.
#'
#' @param assignable An Assignable object, such as a Units or another Clusters object. 
#' @param map A vector that maps units to their cluster ids. The number of
#'   cluster ids must match the number of units in the design.
#' @param covariates An optional data.frame object that contains cluster level covariates.
#' @return A Clusters object.

clusters <- function(assignable, map, covariates = NULL) {
  cs <- nclusters(map) 
  if (number_of_assignments(assignable) != cs) {
    stop("Map vector and assignable objects don't agree on the number of clustered units.")
  }
  new("Clusters", assignable, map = as.integer(map), covariates = covariates)
}

setMethod("number_of_units", signature = c("Clusters"), function(x) { number_of_units(x@.Data) })
setMethod("number_of_assignments", signature = c("Clusters"), function(x) { nclusters(x@map) })
