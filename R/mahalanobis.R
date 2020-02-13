################################################################################
## Functions related to calculating a Mahalanobis distance, moments, etc.
################################################################################

## For each variable in `x` computed the appropriately weighted sum.
##
## @param design A randomized design object
## @param x A covariate matrix (n by k)
## @param z A treatment assignment vector of the length n
## @return A vector of length k of strata weighted sums of each variable
manifest_variable_sums <- function(design, x, z) { UseMethod("manifest_variable_sums") }

## Compute the covariance matrix of the sums returned in `manifest_variable_sums`
##
## @param design A randomized design object
## @param x A covariance matrix (n by k)
## @return A k by k matrix of covariances
manifest_variable_covariance <- function(design, x) { UseMethod("manifest_variable_covariance")}


## Virtual class for MahalanobisDistance objects
setClass("MahalanobisDistance", contains = "numeric")

## First order approximations of Mahalanobis statistics
##
## @slot EM The expected value of the statistic
setClass("ChisquareApproximation",
         contains = "MahalanobisDistance",
         slots = c(EM = "numeric"))

## @slot VarM The variance of the Mahalanobis statistic
## @slot CovT2 Provides the covariance terms Cov(T_i^2 T_j^2) where M = \sum T_k^2
setClass("SecondOrderChisquareApproximation",
         contains = "ChisquareApproximation",
         slots = c(EM = "numeric",
                   VarM = "numeric",
                   CovT2 = "matrix"))

## Compute Mahalanobis distance statistic for two groups
##
## @param x An object giving the design and covariates
## @param z A treatment assignment vector
## @return A MahalanobisDistance object contains both the distance and components useful for calculating moments and other approximations
mahalanobis_distance <- function(x, z) { UseMethod("mahalanobis_distance") }

## If `x` is a matrix, we assume that that is a complete random assignment design
mahalanobis_distance.matrix <- function(x, z) {
    n <- nrow(x)
    d <- create_stratified_design(factor(rep(1, n)), z = z)
    mahalanobis_distance(d, z)
}

## If we can, we see if we can rotate the covariates and then compute mahalanobis_distance
mahalanobis_distance.default <- function(x, z) {
    mahalanobis_distance(rotate_covariates(x), z)
}

mahalanobis_distance.DesignRotatedCovariates <- function(x, z) {
    J <- toJ(x@design, z)
    JTx <- t(J) %*% x
    new("ChisquareApproximation", sum(JTx^2),
        EM  = ncol(x@Rotation))
}
