### Tests we don't want to have run on CRAN   ###
### (ordinarily because the test will throw   ###
### a warning or error during CRAN-checking)  ###
context("Not for CRAN")

###    tests relocated from test.utils.R      ###

## context("Utilities, external dependencies etc")

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

## context("Pseudoinversion via XtX_pseudoinv_sqrt()")

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


## Our plotfun that aims to use RSVGTipsDevice
plotfun_RSVGTD <- function(x, segments, shapes, colors, segments.args, points.args, offset, tiptext) {
  n <- dim(x)[1]
  nstrat <- dim(x)[2]
  ypos <- n:1 + offset

  tts <- "devSVG" == names(dev.cur())[1] && requireNamespace("RSVGTipsDevice")

  if (segments && dim(x)[2] > 1) {
    bnds <- t(apply(x, 1, range))
    do.call(graphics::segments,
            append(list(x0 = bnds[,1],
                        y0 = ypos,
                        x1 = bnds[,2],
                        y1 = ypos),
                   segments.args))
  }

  for(i in 1:nstrat) {

    for (j in seq_along(ypos)) {

      if (tts) {
        # note that these indices are reversed versus convention [i, j, k] notation
        # i is strata (the columns of our tiptext object)
        # j is the variable (the rows of the tips)
        if (dim(tiptext)[3] == 2) {
          RSVGTipsDevice::setSVGShapeToolTip(tiptext[j, i, 1], tiptext[j, i, 2])
        }
        if (dim(tiptext)[3] == 1) {
          RSVGTipsDevice::setSVGShapeToolTip(tiptext[j, i, 1])
        }
      }

      do.call(graphics::points,
              append(list(x[j, i],
                          ypos[j],
                          pch = shapes[j, i],
                          col = colors[j, i]),
                     points.args))
    }
  }


  axis(2, labels = rownames(x), at = ypos, las = 2, tick = FALSE)

  return(offset + n + 1)
}

test_that("Plotting using RSVGTips", {

  # this is just an existence proof: it should go through without errors
  set.seed(20140137)
  
  if(.Platform[['OS.type']]!='windows' && # #71: As of now there are errors in Windows build of
     requireNamespace('RSVGTipsDevice', quietly= TRUE) #the RSVGTipsDevice package; only consider running this on other platforms
     ) {
    require("RSVGTipsDevice")
    f <- tempfile()

    x <- data.frame(z  = rep(c(TRUE, FALSE), 50),
                    x1 = rnorm(100),
                    x2 = sample(c("A", "B", "C"), 100, replace = T),
                    x3 = sample(c("X", "Y", "Z"), 100, replace = T),
                    x4 = sample(c(T,F), 100, replace = T),
                    x5 = sample(c("A", "B", "C"), 100, replace = T),
                    x6 = sample(c("X", "Y", "Z", "W"), 100, replace = T))

    xb <- balanceTest(z ~ x1 * x2 * x3 + strata(x4) + strata(x5), data = x, report = 'all')
    xb$results[, "std.diff", 2] <- xb$results[, "std.diff", 2] * 2

    devSVGTips(paste0(f, "1.svg"), height = 8, width = 8)

    plot(xb, plotfun=plotfun_RSVGTD)

    dev.off()


    xb2 <- balanceTest(z ~ x1 * x2 * x3 + strata(x5), data = x, report = 'all')

    devSVGTips(paste0(f, "2.svg"), height = 8, width = 8)

    plot(xb2, plotfun=plotfun_RSVGTD)

    dev.off()

    xb3 <- balanceTest(z ~ x1 * x2 * x3 + strata(x4) + strata(x5), data = x, report = 'all')
    xb3$results[, "std.diff", 1] <- xb$results[, "std.diff", 1] * 2
    xb3$results[, "std.diff", 2] <- xb$results[, "std.diff", 2] * 4

    devSVGTips(paste0(f, "3.svg"), height = 8, width = 8)
    plot(xb3, plotfun=plotfun_RSVGTD)
    dev.off()
  }
  expect_true(TRUE)

})

test_that("HB08 and HB08_2016 flag degenerate statistics", {

  set.seed(0303022134)

  # pairs
  s <- 20
  n <- 2 * s
  z <- rep(c(0,1), s)
  b <- rep(1:s, each = 2)

  ## theory indicates that the statistic will be degenerate when r - (n - s) = 0
  ## where r is the rank of the matrix in the quadratic form of the test statistic

  x <- replicate(n - s, runif(n, 0, 100))

  colnames(x) <- paste0("x", 1:(n-s))


  df_good <- data.frame(z = z, x = x[, -1], b = b, '(weights)' = 1, check.names = FALSE)
  designs_good <- RItools:::makeDesigns(z ~ . + strata(b) - 1, data = df_good)
  designs_good <- as(designs_good, "StratumWeightedDesignOptions")
  designs_good@Sweights <- RItools:::DesignWeights(designs_good)
  aligned_good <- RItools:::alignDesignsByStrata("b", designs_good)

  expect_silent(HB08(aligned_good))
  expect_silent(HB08_2016(aligned_good))
  
  df_bad <- data.frame(z = z, x = x, b = b, '(weights)' = 1, check.names = FALSE)
  designs_bad <- RItools:::makeDesigns(z ~ . + strata(b) - 1, data = df_bad)
  designs_bad <- as(designs_bad, "StratumWeightedDesignOptions")
  designs_bad@Sweights <- RItools:::DesignWeights(designs_bad)
  aligned_bad <- RItools:::alignDesignsByStrata("b", designs_bad)

  expect_warning(HB08(aligned_bad), "degenerate")
  expect_warning(HB08_2016(aligned_bad), "degenerate")
  
  expect_warning(balanceTest(z ~ . + strata(b), data = df_bad, inferentials.calculator = RItools:::HB08), "degenerate")
  expect_warning(balanceTest(z ~ . + strata(b), data = df_bad, inferentials.calculator = RItools:::HB08_2016), "degenerate")

  skip_on_cran()
  ## `.` includes b as a covariate, and b is constant within its own strata, so
  ## the numerical-stability screen now (correctly) drops it with a message.  The
  ## point of this test is that a NON-degenerate design raises no DEGENERATE
  ## warning, so assert the absence of a warning (the screen message is expected).
  expect_warning(suppressMessages(balanceTest(z ~ . + strata(b), data = df_good, inferentials.calculator = RItools:::HB08)), regexp = NA)
  expect_warning(suppressMessages(balanceTest(z ~ . + strata(b), data = df_good, inferentials.calculator = RItools:::HB08_2016)), regexp = NA)


})

### ACAT omnibus end-to-end null calibration (Monte Carlo; slow) ###
## test.cauchycomb.R checks that acat_pvalue() is uniform under the null at the
## p-value level.  This complements it one level up: the cauchy_comb_p column
## balanceTest() reports must itself be ~Uniform(0,1) under the within-stratum
## null, since the per-covariate p-values feeding it are.  That validity is what
## lets us report ACAT as the omnibus in the degenerate regime; it should already
## hold (GREEN).
test_that("balanceTest cauchy_comb_p is ~uniform under the within-stratum null", {
  skip_on_cran()
  set.seed(2026)
  S <- 25
  ps <- replicate(600, {
    loc <- matrix(rnorm(S * 3), S, 3)
    dif <- matrix(rnorm(S * 3), S, 3)               # no systematic tilt: null true
    X <- matrix(0, 2 * S, 3)
    X[seq(1, 2 * S, 2), ] <- loc + dif / 2          # treated rows
    X[seq(2, 2 * S, 2), ] <- loc - dif / 2          # control rows
    d <- data.frame(z = rep(c(1L, 0L), S), s = factor(rep(seq_len(S), each = 2)),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])
    balanceTest(z ~ x1 + x2 + x3 + strata(s), data = d,
                cauchy.combination = TRUE)$overall["s", "cauchy_comb_p"]
  })
  expect_equal(mean(ps < 0.05), 0.05, tolerance = 0.03)
  expect_equal(mean(ps), 0.5, tolerance = 0.05)
})

### Tests to write...
##test_that("alignDesigns properly tracks UnitWeights vs NotMissing",{})
##test_that("",{})
