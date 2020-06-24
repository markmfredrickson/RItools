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
    (toZ(z) - pZ) 
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
    V <- as.matrix(s %*% diag((p * (sn1 - 1) / (sn - 1) - p^2), ncol = k, nrow = k)  %*% t(s))

    diag(V) <-  as.vector(s %*% ((1 - p) * p))

    ## Since t(x) %*% V %*% x is symmetric, it is diagonalizable as Q D Q^T
    ## with inverse Q^T D^{-} Q (with - indicating any 0 entries are still zero, otherwise 1/d)
    ## the matrix square root of the inverse is therefore Q^T D^{-1/2}
    ## we have a utility for that computation (where the "X" here is x %*% V^{1/2}):
    ## Note: while V is sparse, the result should be dense so there is no harm in the cast to a dense type
    ## and (as of this note) there is not an SVD method for the sparse matrix.
    ## Note on the tolerance arg: in an empirical example, I noticed tha the default in the function
    ## set too many to zero because the first PC had an absolutely giant proportion of the variance
    rotation <- XtX_pseudoinv_sqrt(as.matrix(t(x) %*% V %*% x), mat.is.XtX = TRUE, tol = sqrt(.Machine$double.eps))

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
    euclidean_squared_covariance(design, rotate_covariates(design, covariates))
}

## Method for stratified designs
euclidean_squared_covariance.StratifiedDesign <- function(design, covariates) {
    first_order <- strata_t_covariance_matrices(design, covariates)
    second_order <- strata_t2_covariance_matrices(design, covariates)

    s <- dim(first_order)[1]
    k <- dim(first_order)[2]

    if (s == 1) {
        return(strata_matrix_sum(second_order))
    }

    ## at this point, s >= 2
    tmp <- matrix(0, k, k)
    for (i in 1:(s-1)) {
        for (j in (i+1):s) {
            tmp <- tmp + first_order[i,,] * first_order[j,,]
        }
    }

    return(strata_matrix_sum(second_order) + 4 * tmp)
}

mahalanobis_distribution.StratifiedDesign <- function(x, rotated) {
    covariance_matrix <- as.matrix(euclidean_squared_covariance(x, rotated))

    new("SecondOrderChisquareApproximation",
        EM = ncol(rotated),
        VarM = sum(covariance_matrix),
        CovT2 = covariance_matrix )
}

## Helper to preprocess strata cy by centering and scaling.
## Name comes from the fact we usually call such matrices \tilde X
tilde_maker <- function(design, covariates) {

    ## we will be multiplying by J = Z - P(Z = 1)
    ## knowing that P(Z = 1) = n1/n (by stratum), rearranging terms
    ## shows that T = x'J = (x - \bar x ) ' Z, again all by strata

    strata_means <- as.matrix(t(design@Units) %*% covariates) / design@Count
    centered <- covariates - design@Units %*% strata_means


    return(as.matrix(centered))
}

## Compute the first order covariance matrices
strata_t_covariance_matrices <- function(design, covariates) {
    xtilde <- tilde_maker(design, covariates)

    mu11 <- strata_pairwise_means(design, xtilde)

    ## the quantity $n(N - n)$ shows up frequently
    N <- design@Count
    n1n0 <- design@Treated * (N - design@Treated)
    coef1 <- n1n0 / (N - 1)

    if (length(design@Count) == 1) {
        v <- dim(xtilde)[2]
        mu11 <- array(mu11, c(1, v, v))
    }


    coef1 * mu11
}

## Helper to get second order covariance matrices
strata_t2_covariance_matrices <- function(design, covariates) {

    xtilde <- tilde_maker(design, covariates)

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
        mu11 <- array(mu11, c(1, v, v))
        mu22 <- array(mu22, c(1, v, v))
        mu2_mu2 <- array(mu2_mu2, c(1, v, v))
    }

    ## the quantity $n(N - n)$ shows up frequently
    N <- design@Count
    n1n0 <- design@Treated * (N - design@Treated)

    ## For any estimate with the same order, we need the same coefficients when
    ## computing the strata level expected values

    ## for strata with fewer than 4 units, we need to special case the
    ## coefficients so that result in the multiplication is equal to mu22
    coef1 <- n1n0 / (N - 1)
    coef2 <- ifelse(N > 3, coef1 / ((N - 2) * (N - 3)), 1)
    coef2a <- ifelse(N > 3, (N * (N + 1) - 6 * n1n0), 1)
    coef2b <- ifelse(N > 3, N * (N - 1 - n1n0), 0)

    ## strata expected values E(T_k^2 T_j^2)
    mean_22 <- coef2 * (coef2a * mu22 - coef2b * (2 * mu11^2 + mu2_mu2))

    ## strata products of expected values E(T_k^2) E(T_j^)
    mean2_mean2 <- coef1^2 * mu2_mu2

    ## now compute the (per stratum) covariance matrices $E(T^2 [T^2]') - E(T^2) E(T^2)'
    strata_covariance_array <- mean_22 - mean2_mean2

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
    s <- dim(a)[1]
    tmp <- a[1,,]
    if (s > 1) {
        for (j in 2:s) {
            tmp <- a[j,,] + tmp
        }
    }

    return(tmp)
}

