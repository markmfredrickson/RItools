### Tests we don't want to have run on CRAN   ###
### (ordinarily because the test will throw   ###
### a warning or error during CRAN-checking)  ###

###    tests relocated from test.utils.R      ###

context("Utilities, external dependencies etc")

test_that("data.table options issue #69", {

  if (suppressMessages(suppressWarnings(require(data.table)))) {
    data(nuclearplants)
    f <- function() 1
    expect_equal(withOptions(list(), f), 1)
  }

})

test_that("survival::strata() still conforms to our expectations",
{
    facA <- rep(c("a", "b", "c"), each=2)
    facB <- rep(c("A", "B", "C"), 2)
    facC <- survival::strata(facA, facB)
    expect_equal(6, nlevels(facC))
})

context("Pseudoinversion via XtX_pseudoinv_sqrt()")

test_that("proper inversion in full rank case",{
    basis <- cbind(1/sqrt(3), poly(rnorm(3), degree=2))
    expect_gt(abs(det(basis)),.Machine$double.eps) # full rank
    mat <- t(basis) %*% diag((2:0)) %*% basis
    m2 <- XtX_pseudoinv_sqrt(mat)
    expect_equivalent(basis %*% tcrossprod(m2) %*% t(basis), diag(c(1/2^2, 1, 0)) )
})


##load("tricky_rectangular_matrix.rda")
##Make a matrix that is singular
##tricky_rectangular_matrix <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
##                                      0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1), 9, 4)
n <- 14 
tricky_rectangular_matrix <- matrix(NA, n, n) 
for (i in 1:n) for (j in 1:n) tricky_rectangular_matrix[i,j] <- 1/(i+j-1) 

test_that("answers match MASS::ginv() under near rank deficiency",{
    XtX  <- crossprod(tricky_rectangular_matrix)
    pinv_XtX  <-  MASS::ginv(XtX)
    pinv_sqrt_XtX  <- XtX_pseudoinv_sqrt(tricky_rectangular_matrix,
                                         tol=.Machine$double.eps^0.25)
    expect_equivalent(pinv_XtX, tcrossprod(pinv_sqrt_XtX))
})
