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
    (toZ(z) - pZ) / pZ 
}


## Method for stratified designs
rotate_covariates.StratifiedDesign <- function(design, x) {
    ## convenient shorthand
    s <- design@Units
    sn <- design@Count
    sn1 <- design@Treated
    sn0 <- sn - sn1
    k <- length(sn)
    n <- sum(sn)

    ## strata level probabilities
    p <- sn1 / sn # pi

    ## TODO: I would like to write this as a smaller k by k matrix, where k is number of strata 
    ## but I couldn't quite figure out how to do that. It might not be possible.

    ## the matrix E(JJ'), J_i = (Z_i - pi_i) / pi_i = (Z_i - n1/n) / (n1 / n)
    V <- s %*% diag((p * (sn1 - 1) / (sn - 1) - p^2) / p^2, ncol = k, nrow = k)  %*% t(s)

    ## For some reason, this doesn't work in some cases
    ## diag(V) <- 1 / as.vector(s %*% p_one_p)
    diag_V <-  as.vector(s %*% ((1 - p) / p))
    for (i in 1:n) {
        V[i,i] <- diag_V[i]
    }

    ## Since t(x) %*% V %*% x is symmetric, it is diagonalizable as Q D Q^T
    ## with inverse Q^T D^{-} Q (with - indicating any 0 entries are still zero, otherwise 1/d)
    ## the matrix square root of the inverse is therefore Q^T D^{-1/2}
    ## we have a utility for that computation (where the "X" here is x %*% V^{1/2}):
    ## Note: while V is sparse, the result should be dense so there is no harm in the cast to a dense type
    ## and (as of this note) there is not an SVD method for the sparse matrix.
    rotation <- XtX_pseudoinv_sqrt(as.matrix(t(x) %*% V %*% x), mat.is.XtX = TRUE)

    new("DesignRotatedCovariates", x %*% rotation,
        Design = design,
        Rotation = rotation)
}

## Method for stratified designs
manifest_variable_sums.StratifiedDesign <- function(design, x, z) {

    n1_over_n <- as.vector(design@Units %*% (design@Treated / design@Count))

    ssn <- sparseToVec(t(matrix(toZ(z), ncol = 1) - n1_over_n) %*% x, column = FALSE)
    names(ssn) <- colnames(x)

    return(ssn)
}

## Method for stratified designs
manifest_variable_covariance.StratifiedDesign <- function(design, x) {
    ## TODO: is the returned thing just the rotation from rotate_covariates?

    ## the fact that this is duplicated suggests it should be a function or pre-computed in an object
    p_ <- ncol(x)
    s_ <- ncol(design@Units)

    ## this next block first creates the n_ * p_^2 matrix
    ## of 2nd-order monomials in columns of x-tilde, then
    ## immediately sums each of these within each of s_ strata,
    ## resulting in a n_ * p_^2 matrix.
    xt_df  <- as.data.frame(x) # to get rep() to treat as a list
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

## Method for stratified designs
t_squared_covariance.StratifiedDesign <- function(design, covariates) {
    strata_mats <- strata_covariance_matrices(design, covariates)
    strata_matrix_sum(strata_mats)
}

## Helper to get per strata covariance matrices
strata_covariance_matrices <- function(design, covariates) {
    ## overall we are computing (for each strata)
    ## E(T^2 T^2') - E(T^2) E(T^2)'
    ## where T = (T_1, ..., T_k)'

    rotated <- rotate_covariates(design, covariates)

    ## rotated assumes we will be multiplying by J = (Z - P(Z = 1)) / P(Z = 1) 
    ## knowing that P(Z = 1) = n1/n (by stratum), rearranging terms
    ## shows that T = x'J = [(x - \bar x ) / (n1/n)] ' Z, again all by strata

    strata_means <- as.matrix(t(design@Units) %*% rotated) / design@Count

    rotated_centered <- rotated - design@Units %*% strata_means

    ## now multiply through by n/n1
    xtilde <- rotated_centered * as.vector(design@Units %*% (design@Count / design@Treated))
    xtilde <- as.matrix(xtilde) ## make sure this is dense

    ## Finucan uses a subscript notation for, eg., \mu_{ab} = N^{-1} \sum x_{ia} x_{ib} 
    ## with \sigma_a^2 = \mu_{aa}, etc

    mu11 <- strata_pairwise_means(design, xtilde)
    mu22 <- strata_pairwise_means(design, xtilde^2)

    ## the last piece we need is the products E(T_k^2) E(T_j^2) per strata
    mu2 <- as.matrix(t(design@Units) %*% xtilde^2) / design@Count
    mu2_mu2 <- pairwise_products(mu2) 

    ## handling single stratum case:
    if (length(design@Count) == 1) {
        v <- dim(xtilde)[2]
        mu11 <- array(mu11, c(v, v, 1))
        mu22 <- array(mu22, c(v, v, 1))
        mu2_mu2 <- array(mu2_mu2, c(v, v, 1))
    }

    ## the quantity $n(N - n)$ shows up frequently
    N <- design@Count
    n1n0 <- design@Treated * (N - design@Treated)

    ## For any estimate with the same order, we need the same coefficients when
    ## computing the strata level expected values
    ## TODO: fix of strata with fewer then 2 or 4 units (respectively)
    coef1 <- n1n0 / (N - 1)
    coef2 <- coef1 / ((N - 2) * (N - 3))
    coef2a <- (N * (N + 1) - 6 * n1n0)
    coef2b <- N * (N - 1 - n1n0)

    ## now we strata expected values
    mean_22 <- coef2 * (coef2a * mu22 - coef2b * (2 * mu11^2 + mu2_mu2))

    ## now compute the (per stratum) covariance matrices $E(T^2 [T^2]') - E(T^2) E(T^2)'
    ## E(T^2) = 1, so
    strata_covariance_array <- mean_22 - 1

    return(strata_covariance_array)
}

## @param x A matrix of n by k
## @return An array of n by k by k, with each [i, j, l] entry being x[i, j] * x[i, l]
pairwise_products <- function(x) {
    n <- dim(x)[1]
    k <- dim(x)[2]

    ## the apply will get all the column wise products, but the ordering isn't quite what we want
    aperm(perm = c(1, 3, 2),
          array(apply(x, 2, function(xx) { as.vector(x * xx)}),
                dim = c(n, k, k)))
}

## Returns strata level means of \sum_i x_{ij} x_{ik}
## @param d A Stratified Design object
## @param x An array of n by k
## @return An array of s by k by k indicating strata level means of x_i * x_j
strata_pairwise_means <- function(d, x) {
    xx <- pairwise_products(x)
    sx <- apply(xx, 2:3, function(xi) { t(d@Units) %*% xi / d@Count })
    return(sx)
}

## Sum matrices over strata
## @param a s by j by k matrix
## @return A j by k matrix produced by summing the other component matrices
strata_matrix_sum <- function(a) {
    s <- dim(a)[3]
    tmp <- a[,,1]
    if (s > 1) {
        for (j in 2:s) {
            tmp <- a[,,j] + tmp
        }
    }

    return(tmp)
}
