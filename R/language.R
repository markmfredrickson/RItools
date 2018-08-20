################################################################################
## A language of designs
##
## Specifying designs is a key aspect of analyzing randomized trials. These
## objects encapsulate critical design elements and then reduce those to
## appropriate quantities during computation steps.
##
################################################################################


## Generic functions that should have methods implemented by Design objects
## TODO: investigate "group generic" functions -- can these be defined by
## programmers or are they only for internal implementations

#' The number of treatment levels for a Design object.
#'
#' @param x The Design object.
#' @return The number of different treatment levels in the design.
setGeneric("treatment_levels", def = function(x) { NA })


#' Update subgrouping information for the design.
#'
#' @param x The Design object.
#' @param subgroup A factor of length n indicating which subgroup each unit is
#'   in. NA values indicate the unit is not a member of any subgroup to be
#'   analyzed (e.g., is missing on an outcome).
#' @return A new design with the subgrouping updated.
setGeneric("subgroup", def = function(design, subgroup) { stop("Not yet implemented.") })

## a little helper to grab different subgroups and do something to them
subgroup_helper <- function(design, f) {
  ll <- lapply(levels(design@subgroup), function(lvl) {
    idx <- design@subgroup == lvl
    idx[is.na(idx)] <- FALSE
    f(idx)
  })
  names(ll) <- levels(design@subgroup)
  return(ll)
}

## TODO: should this have the HT name, or something else?
## Investigate glht to see what they call them. It is basically that kind of object.
setClass("HorvitzThompsonTotals",
         slots = list(levels = "factor", subgroups = "factor", estimates = 'numeric', covariance = 'matrix'))

## default implementation for a generic. keeping it as separate function for
## debugging/testing
.horvitz_thompson <- function(design, y, z) {
    stop("Not yet implemented")
}

#' Horvitz-Thompson estimator of a total
#'
#' A generic function for generating Horvitz-Thompson inverse probability
#' estimates of finite population totals. This function provides a generic
#' solution for any "Design" object using first and second order probabilities
#' of assignment. Particular designs can provide more efficient or numerically
#' stable implementations if needed.
#'
#' @param design A "Design" object that pertains to n units.
#' @param x The observed variable to estimate (length n)
#' @param z An integer vector of length n indicating to which condition each unit is assigned.
#'
#' @return A vector of estimates of the totals for all levels in the study and an estimate of the covariance matrix for the totals. TODO: investigate the GLHT package to see what kind of objects they take. 

setGeneric("horvitz_thompson", def = function(design, y, z) { stop("Not yet implemented.")})

# TODO: might want to validate input here for y and z for all functions (length, levels, etc)

## Virtual classes to for the designs, which have units and clusters, but organize them in some way.
## The AbstractRCT class just show how things are randomized and can be used to generate randomizations.
## The ContcreteRCT takes an abstract RCT and combines it with a randomization factor and something 
setClass("AbstractRCT")
setClass("ConcreteRCT", contains = "AbstractRCT", slots = c(target = "Assignable", assignment = "factor"))

makeConcrete <- function(abstractRCT, target, assignment) {
    nassignable <- number_o
    new("ConcreteRCT", abstractRCT, )
}

setClass("AbstractCompleteRandomAssignment",
         contains = "AbstractRCT",
         slots = list(levels = "integer"))

#' Create design object for complete random assignment
#'
#' For J categories, assign a fixed number of clusters to each category.
#'
#'
#' @param treatments A vector of length J, indicating how many clusters will be
#'   assigned to each treatment condition. The sum of all levels should be equal
#'   to the total number of clusters (n).
#'
#' @return A "CompleteRandomAssignment" object

complete_random_assignment <- function(treatments) {
  new("CompleteRandomAssignment", levels = as.integer(treatments), subgroup = as.factor(rep(1, sum(treatments))))
}

setMethod("treatment_levels", signature = "CompleteRandomAssignment", function(x) { length(x@levels) })
setMethod("number_of_units", signature = "CompleteRandomAssignment", function(x) {
  tmp <- as.vector(table(x@subgroup, useNA = 'ifany'))
  names(tmp) <- levels(x@subgroup)
  return(tmp)
})


#' @import magic adiag
.horvtiz_thompson.CompleteRandomAssignment <- function(design, y, z) {
  ## for convenience of notation
  n_per_group <- number_of_units(design)
  n <- sum(n_per_group)
  z <- as.factor(z)

  J <- treatment_levels(design)
  nk <- design@levels
  G <- nlevels(design@subgroup)

  if (length(y) != n) {
    stop("Length of outcome not equal to the number of units.")
  }

  if (length(z) != n) {
    stop("Length of the treatment assignment vector not equal to the number of units.")
  }

  f <- function(idx) {
    yy <- y
    yy[!idx] <- 0

    ## TODO: should we try for a sparse model matrix? or do this with loops in C?
    ## another alternative would be to use tapply and mean
    membership <- model.matrix(~ z - 1)

    ytot <- as.vector(t(membership) %*% yy)
    ybar <- ytot / nk

    estimates <- n * ybar

    ## again, might be able to optimize this
    y2tot <- as.vector(t(membership) %*% yy^2)
    s2 <- (y2tot -  nk * ybar^2) / (nk - 1)
    vars <- n * (n - nk) / nk * s2

    sigmahats <- (n - 1) / n * s2

    covs <- -(0.5) * n^2/(n - 1) * outer(sigmahats, sigmahats, `+`)
    diag(covs) <- vars

    list(levels = colnames(membership), estimates = estimates, covariance = covs)
  }

  perSubgroup <- subgroup_helper(design, f)

  new("HorvitzThompsonTotals",
      levels = as.factor(sapply(perSubgroup, function(h) { h$levels })),
      subgroups = as.factor(sapply(names(perSubgroup), function(nm) { rep(nm, J) })),
      estimates = unlist(lapply(perSubgroup, function(h) { h$estimates })),
      covariance = Reduce(function(x, y) { magic::adiag(x, y, pad = NA) },
                          lapply(perSubgroup, function(l) { l$covariance })))
}

setMethod("horvitz_thompson", signature = c("CompleteRandomAssignment", "numeric", "numeric"), .horvtiz_thompson.CompleteRandomAssignment)

setClass("BlockedRandomAssignment",
         contains = "Assignable", 
         slots = list(blocks = "matrix", map = "integer"))

#' Complete random assignment within blocks.
#'
#' For M blocks, each assigns some number of clusters to J categories.
#'
#' @param blocks A J by M matrix of the number of units within each block
#'   (column) assigned to each treatment condition (row).
#' @param map A vector of block assignments that will map any vector of unit data to block membership.
#'
#' @return TODO: get a template for these functions to avoid duplicating docs. Or use a common help page.

blocked_random_assignment <- function(blocks, map) {
  storage.mode(blocks) <- "integer"
  if (identical(rowSums(blocks), table(map))) {
   stop("Map vector and blocks matrix do not agree on the number of units per block.")
  }
  new("BlockedRandomAssignment", blocks = blocks, map = as.integer(map), subgroup = as.factor(rep(1, sum(blocks))))
}

setMethod(number_of_units, signature = c("BlockedRandomAssignment"), function(x) { as.vector(table(x@subgroup)) })
setMethod(treatment_levels, signature = c("BlockedRandomAssignment"), function(x) { ncol(x@blocks) })


.horvtiz_thompson.BlockedRandomAssignment <- function(design, y, z)  {
  ## for convenience of notation
  n <- number_of_units(design)
  J <- treatment_levels(design)
  blks <- design@blocks
  B <- nrow(blks)

  estimates <- numeric(J)
  vcv <- matrix(0, ncol = J, nrow = J)

  ht <- NULL
  for (b in 1:B) {
    ht <- horvitz_thompson(complete_random_assignment(blks[b, ]), y[design@map == b], z[design@map == b])
    estimates <- estimates + ht@estimates
    vcv <- vcv + ht@covariance
  }

  ## use the last value of ht to get the ordering of levels and subgroups,
  ## as they should be the same across the blocks
  new("HorvitzThompsonTotals", levels = ht@levels, subgroups = ht@subgroups, estimates = estimates, covariance = vcv)
}

setMethod("horvitz_thompson", signature = c("BlockedRandomAssignment", "numeric", "numeric"),
          .horvtiz_thompson.BlockedRandomAssignment)


## Helper class for things like clusters and subgroups that wrap existing designs

setClass("Clusters",
         contains = "Assignable", 
         slots = list(map = "factor"))

#' Clustered treatment assignment mechanisms.
#'
#' Units are arranged into clusters, which form the basis of randomization.
#' Units are mapped to clusters when computing quantities. The clusters require
#' their own design object.
#'
#' @param design A Design object (e.g., a CompleteRandomAssignment) that uses
#'   the *clusters* as the randomization unit.
#' @param map A vector that maps units to their cluster ids. The number of
#'   cluster ids must match the number of units in the design.
#' @return A ClusteredAssignment object.

clustered_assignment <- function(design, map) {
  cs <- length(unique(map))
  if (number_of_units(design) != cs) {
    stop("Map vector and cluster level design object do not agree on the number of clusters.")
  }
  new("ClusteredAssignment", design = design, map = as.integer(map), subgroup = as.factor(rep(1, length(map))))
}

setMethod(number_of_units, signature = c("ClusteredAssignment"), function(x) { as.vector(table(x@subgroup)) })

setMethod("horvitz_thompson",
           signature = c("ClusteredAssignment", "numeric", "numeric"),
           function(design, y, z) {
             clusterZ <- table(design@map, z) > 0

             if (any(rowSums(clusterZ) != 1)) {
                 stop(paste("Some clusters contain multiple treatment assignments (",
                            paste(which(rowSums(clusterZ) != 0), collapse = ","),
                            ")",
                            sep = ""))
             }
             zz <- apply(clusterZ, 1, which.max)

             f <- function(idx) {
                 yy <- y
                 yy[!idx] <- 0

                 ysum <- as.numeric(tapply(yy, design@map, sum))

                 return(horvitz_thompson(design@design, ysum, zz))
             }

             sgps <- subgroup_helper(design, f)

             ## TODO: factor this code out of multiple HT functions
             new("HorvitzThompsonTotals",
                 levels = as.factor(unlist(lapply(sgps, function(x) { x@levels }))),
                 subgroups = as.factor(rep(names(sgps), each = length(sgps[[1]]@estimates))),
                 estimates = unlist(lapply(sgps, function(h) { h@estimates })),
                 covariance = Reduce(function(x, y) { magic::adiag(x, y, pad = NA) },
                                     lapply(sgps, function(l) { l@covariance })))


           })

subgroup <- function(design, subgroup) {
  if (length(subgroup) != sum(number_of_units(design))) {
    stop("Subgroup indicator does not agree with the number of units in the design.")
  }

  design@subgroup <- factor(interaction(design@subgroup, subgroup))
  return(design)
}
