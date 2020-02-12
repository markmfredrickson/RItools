## Stratified design of x units into k strata with fixed numbers of units and
## treated per strata
##
## @slot Unit A n by k (sparse) matrix with a single 1 in each row. Rows
##   indicated the randomized units (clusters in cluster randomized trials).
##   Columns indicate stratum membership.
## @slot Counts A k vector of units in each stratum
## @slot Treated A k vector of the number treated in each stratum.
## @slot Weights A k vector of strata weights.
## @slot JCoefs A n vector of precomputes coefficients to produce J = Z (1 - P(Z = 1)) / (P(Z = 1) P(Z = 0))
setClass("StratifiedDesign", contains = "RandomizedDesign",
         slots = c(
             Units = "matrix.csr",
             Count = "integer",
             Treated = "integer"
         ))

## Create stratified design object.
##
## @param strata A factor of length n, with k levels.
## @param treated (Optional) An integer vector of length k, with names
##   corresponding to levels of strata, indicating treated units. Either `treated`
##   or `z` must be included
## @param z A logical or two level vector of length n indicating if each unit is
##   treated or not.
## @param weights (Optional) A numeric of length n.
## @return An object of type "StratifiedDesign"
create_stratified_design <- function(strata, treated = NULL, z = NULL, weights = NULL) {
    n <- length(strata)
    k <- nlevels(strata)

    units <- SparseMMFromFactor(strata)

    count <- as.vector(table(strata))

    if (is.null(treated) && !is.null(z)) {
        treated <- as.vector(t(units) %*% toZ(z))
    }

    if (!is.null(names(treated))) {
       treated <- treated[levels(strata)]
    }


    new("StratifiedDesign",
        Units   = units,
        Count   = as.integer(count),
        Treated = as.integer(treated))
}

## method for stratified designs
toJ.StratifiedDesign <- function(x, z) {
    pZ <- as.vector(x@Units %*% (x@Treated / x@Count)) ## P(Z = 1)
    (toZ(z) - pZ) / (pZ * (1 - pZ))
}


## Method for stratified designs
rotate_covariates.StratifiedDesign <- function(design, x) {
    ## convenient shorthand
    s <- design@Units
    sn <- design@Count
    sn1 <- design@Treated
    sn0 <- sn - sn1

    ## strata level probabilities
    p <- sn1 / sn # pi
    one_p <- sn0 / sn # 1 - pi
    p_one_p <- p * one_p # pi (1 - pi)

    ## TODO: I would like to write this as a smaller k by k matrix, where k is number of strata 
    ## but I couldn't quite figure out how to do that. It might not be possible.

    ## the matrix E(JJ'), J_i = (Z_i - pi_i) / (pi_i (1 - pi_i)) = (Z_i - n1/n) / (n1 n0 / n^2)
    V <- s %*% diag((p * (sn1 - 1) / (sn - 1) - p^2) / (p^2 * one_p^2))  %*% t(s)
    diag(V) <- 1 / as.vector(s %*% p_one_p)

    ## Since t(x) %*% V %*% x is symmetric, it is diagonalizable as Q D Q^T
    ## with inverse Q^T D^{-} Q (with - indicating any 0 entries are still zero, otherwise 1/d)
    ## the matrix square root of the inverse is therefore Q^T D^{-1/2}
    ## we have a utility for that computation (where the "X" here is x %*% V^{1/2}):
    rotation <- XtX_pseudoinv_sqrt(t(x) %*% V %*% x, mat.is.XtX = TRUE)

    new("DesignRotatedCovariates", x %*% rotation,
        Design = design,
        Rotation = rotation)
}
