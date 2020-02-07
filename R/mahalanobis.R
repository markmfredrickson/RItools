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

manifest_variable_sums.StratifiedDesign <- function(design, x, z) {

    x_tilde <- x * as.vector(design@Units %*% design@Weights)

    n1_over_n <- as.vector(design@Units %*% (design@Treated / design@Count))

    ssn <- sparseToVec(t(matrix(toZ(z), ncol = 1) - n1_over_n) %*% x_tilde, column = FALSE)
    names(ssn) <- colnames(x)

    return(ssn)
}

## Compute the covariance matrix of the sums returned in `manifest_variable_sums`
##
## @param design A randomized design object
## @param x A covariance matrix (n by k)
## @return A k by k matrix of covariances
manifest_variable_covariance <- function(design, x) { UseMethod("manifest_variable_covariance")}

manifest_variable_covariance.StratifiedDesign <- function(design, x) {

    ## the fact that this is duplicated suggests it should be a function or pre-computed in an object
    x_tilde <- x * as.vector(design@Units %*% design@Weights)
    p_ <- ncol(x_tilde)
    s_ <- ncol(design@Units)

    ## this next block first creates the n_ * p_^2 matrix
    ## of 2nd-order monomials in columns of x-tilde, then
    ## immediately sums each of these within each of s_ strata,
    ## resulting in a n_ * p_^2 matrix.
    xt_df  <- as.data.frame(x_tilde) # to get rep() to treat as a list
    xt_covar_stratwise  <- t(design@Units) %*%
        ( as.matrix( as.data.frame(rep(xt_df, each=p_)) ) *
          as.matrix( as.data.frame(rep(xt_df, times=p_)) )
            )
    xt_covar_stratwise  <- as.matrix(xt_covar_stratwise) # s_ * (p_^2)
    xt_covar_stratwise  <- array(xt_covar_stratwise,
                                 dim=c(s_, p_, p_) # s_ * p_ * p_
                                 )
    xt_covar_stratwise  <-
        ifelse(design@Count == 1, 0, (design@Count -1)^(-1) ) *
        xt_covar_stratwise # still s_ * p_ * p_
    ## now we have sample covariances, by stratum.

    ## Combine with factors equal to half the harmonic means of n1 and n0. 
    xt_c_s_scaled  <- xt_covar_stratwise /
        (1/(design@Count - design@Treated) + 1/design@Treated) 
    tcov  <- apply(xt_c_s_scaled, 2:3, sum)

    return(tcov)
}


## Virtual class for MahalanobisDistance objects
setClass("MahalanobisDistance", representation = "numeric")

## First order approximations of Mahalanobis statistics
##
## @slot EM The expected value of the statistic
setClass("ChisquareApproximation",
         representation = "MahalanobisDistance",
         slots = c(EM = "numeric"))

## @slot VarM The variance of the Mahalanobis statistic
## @slot CovT2 Provides the covariance terms Cov(T_i^2 T_j^2) where M = \sum T_k^2
setClass("SecondOrderChisquareApproximation",
         representation = "MahalanobisDistanceChisquareApproximation",
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
    new("ChisquareApproximation", sum(JTx^2)
        EM  = ncol(x@Rotation))
}
