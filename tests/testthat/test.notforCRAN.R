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


context("Diagnose and fix a problem with balanceTest descriptives")
## Evaluate issue with floating point equalities in balanceTest
source("dumpdata.R")
source("moredat.R")
xb0i <- balanceTest(baselineFmlaCluster, data = dat17i, report = "all", p.adjust.method = "none")

test_that("First, that we end  up with an empty matrix when we just wanted to delete one row",{
expect_equal(nrow(descriptives),82)
		  bad <- apply(descriptives[nmvars, group_mean_labs,,drop=FALSE]==1,1,all)
toremove <- match(nmvars[bad], dimnames(descriptives)[["vars"]])
expect_equal(toremove,82)
descriptives2 <- descriptives[-toremove,,,drop=FALSE]
expect_equal(nrow(descriptives2),81)
})


context("Diagnose and fix a problem with balanceTest descriptives")
## Evaluate issue with floating point equalities in balanceTest
test_that("First, that we end  up with an empty matrix when we just wanted to delete one row",{
		  source("moredat.R")
		  baselineFmla <- reformulate(covs3, response = "soldvsnot17")
		  baselineFmlaCluster <- update(baselineFmla, . ~ . + cluster(Q56))
		  xb0i <- balanceTest(baselineFmlaCluster, data = dat17i, report = "all", p.adjust.method = "none")
		  expect_equal(nrow(xb0i$results),0)
		  ## This  next  does not work  even though the objects were exported from with the debug session of the call to balanceTest above.
		  ##Browse[2]> save(nmvars,group_mean_labs,descriptives,origvars,file="objects_from_debug_balanceTest.rda")
		  load("objects_from_debug_balanceTest.rda")
		  ## The problem lines from lines 289--292 in balanceTest.R
		  bad1 <- apply(descriptives[nmvars, group_mean_labs,,drop=FALSE]==1,1,all)
		  toremove1 <- match(nmvars[bad1], dimnames(descriptives)[["vars"]])
		  expect_equal(toremove1,integer(0))
		  descriptives_gone	 <- descriptives[-toremove1,,,drop=FALSE]
		  origvars_gone	 <- origvars[-toremove1]
		  expect_equal(nrow(descriptives_gone),0)
		  expect_equal(dim(origvars_gone),NULL)
		  ## Now showing one inelegant fix
		  bad2 <- apply(descriptives[nmvars, group_mean_labs,,drop=FALSE],1,
			       function(x){ all.equal(x,matrix(rep(1,length(x)),nrow=1),check.attributes=FALSE) })
		  toremove2 <- match(nmvars[bad2], dimnames(descriptives)[["vars"]])
		  expect_equal(toremove2,82)
		  descriptives_ok <- descriptives[-toremove2,,,drop=FALSE]
		  expect_equal(nrow(descriptives_ok),81)
		  origvars_ok <- origvars[-toremove2]
		  expect_equal(origvars_ok,1:81)
		  ## A better fix?
		  groupmeans <- descriptives[nmvars, group_mean_labs,,drop=FALSE]
		  ## sqrt(.Machine$double.eps) is the default tolerance in all.equal()
		  bad3 <- apply(abs(groupmeans - 1)<sqrt(.Machine$double.eps), 1, all)
		  toremove3 <- match(nmvars[bad3], dimnames(descriptives)[["vars"]])
		  expect_equal(toremove3,82)
		  descriptives_ok <- descriptives[-toremove3,,,drop=FALSE]
		  expect_equal(nrow(descriptives_ok),81)
		  origvars_ok <- origvars[-toremove3]
		  expect_equal(origvars_ok,1:81)
})


