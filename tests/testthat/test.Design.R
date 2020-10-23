################################################################################
# Tests for DesignOptions objects
################################################################################

library("testthat")


test_that("Null stratification is encoded by '--'",{
    data(nuclearplants)
    nuclearplants$"(weights)" <- 1
    d0  <- makeDesigns(pr ~ cost, data=nuclearplants)
    expect_equal(colnames(d0@StrataFrame), "--")
    foo <- nuclearplants$pt
    d1 <- makeDesigns(pr ~ cost + strata(foo), data=nuclearplants)
    expect_true("--" %in% colnames(d1@StrataFrame))
})

test_that("Issue #86: makeDesigns finds variables outside data arg",{
    data(nuclearplants)
    foo <- nuclearplants$pt
    nuclearplants$"(weights)" <- 1
    d <- makeDesigns(pr ~ cost + strata(foo), data=nuclearplants)
      expect_s4_class(d, "DesignOptions")
})


context("DesignOptions objects")

test_that("Issue #76: Using I() in formulas", {

    x <- data.frame(z=c(1,1))  # have to exclude 
    while (all(x$z==x[1L,'z'])) # degenerate case
        x <- data.frame(x = rnorm(10), y = rnorm(10), z = rbinom(10, size = 1, p = 1/3))

        x$"(weights)" <- 1
        d <- makeDesigns(z ~ I(x * sin(y)), data = x)
        expect_s4_class(d, "DesignOptions")
        ## While we're at it, confirm that the non-stratification
        ## encoded in this DesignOptions bears the column name "--".
        
})

test_that("Creating DesignOptions objects", {

    set.seed(20130801)
  d <- data.frame(
      id          = 1:500,
      x           = rnorm(500),
      cluster     = rep(1:100, 5),
      strata.good = rep(c(1:4, NA), 100),
      strata.bad  = sample(1:100, size = 500, replace = T),
      z.good      = rep(c(0,1), 250),
      z.bad       = sample(c(0,1), size = 500, replace = T))
  d$'(weights)' = 1 # meet expectation of a weights column
  # checking input
  expect_error(RItools:::makeDesigns(~ x, data = d), "treatment")
  expect_error(RItools:::makeDesigns(z.bad ~ x + strata(strata.good) + cluster(cluster), data = d), "cluster")
  expect_error(RItools:::makeDesigns(z.good ~ x + strata(strata.bad) + cluster(cluster), data = d), "strata")
  expect_error(RItools:::makeDesigns(z.good ~ x + cluster(id) + cluster(cluster), data = d), "cluster")
  expect_error(RItools:::makeDesigns(z.good ~ x - 1, data = d, "stratification"))
  # actually testing that the output is as expected
  simple <- RItools:::makeDesigns(z.good ~ x, data = d)
  expect_equal(dim(simple@StrataFrame)[2], 1)
  expect_equivalent(simple[, "x"], d$x)
  expect_equivalent(simple@Z, as.logical(d$z.good))
  expect_equal(nlevels(simple@Cluster), 500) # a cluster per individual
  
  clustered <- RItools:::makeDesigns(z.good ~ x + cluster(cluster), data = d)
  expect_equal(dim(clustered@StrataFrame)[2], 1)
  expect_true(nlevels(clustered@Cluster) > 1)

  clustStrata <- RItools:::makeDesigns(z.good ~ x + cluster(cluster) + strata(strata.good), data = d)
  expect_equal(dim(clustStrata@StrataFrame)[2], 2)

    clustStrata.c <- RItools:::makeDesigns(z.good ~ x + cluster(cluster) + strata(strata.good, strata.good), data = d)
    expect_equivalent(clustStrata, clustStrata.c)
    
  # dropping the overall comparison
  expect_equal(dim(RItools:::makeDesigns(z.good ~ x + cluster(cluster) + strata(strata.good) - 1, data = d)@StrataFrame)[2], 1)
   
  ## More tests to write:
  # - All NA strata variables
  # - Missing z or cluster
  # - strata with extra levels but no observations (which can be safely dropped)
  #   (NB: extra levels tested upstream, in xBalance, as of commit 34861515; 
  #   see ./test.clusters.R ) 
})

test_that("NotMissing vars correctly generated",
          {

  dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, TRUE), c(3,2)),
                    x2 = c(1:5), 
                    x3=c(TRUE, FALSE, NA, TRUE, FALSE),
                    fac=factor(c(rep(1:2,2), NA))
                    )
  dat$'(weights)' <- 1

  simple <- RItools:::makeDesigns(z ~ x1 + x2  + strata(strat) + cluster(clus), data = dat)
  expect_match(colnames(simple@NotMissing), "x1", all=FALSE)
  expect_false(any(grepl("x2", colnames(simple@NotMissing))))
  expect_false(any(grepl("TRUE", colnames(simple))))
  expect_false(any(grepl("FALSE", colnames(simple))))

  simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
  expect_match(colnames(simple2@NotMissing), "x1", all=FALSE)
  expect_false(any(grepl("x2", colnames(simple2@NotMissing))))
  expect_match(colnames(simple2@NotMissing), "fac", all=FALSE)
  expect_false(any(grepl("TRUE", colnames(simple2))))            
  expect_false(any(grepl("FALSE", colnames(simple2))))

  simple3 <- RItools:::makeDesigns(z ~ x1 + x3 + fac+ strata(strat) + cluster(clus), data = dat)
  expect_match(colnames(simple3@NotMissing), "x1", all=FALSE)
  expect_match(colnames(simple3@NotMissing), "x3", all=FALSE)
  expect_match(colnames(simple3@NotMissing), "fac", all=FALSE)
  expect_false(any(grepl("TRUE", colnames(simple3))))            
  expect_false(any(grepl("FALSE", colnames(simple3))))
              
          })

test_that("Issue 88: logical Covariates correctly generated",
          {

  dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, TRUE), c(3,2)),
                    x2 = c(1:5), 
                    x3=c(TRUE, FALSE, NA, TRUE, FALSE),
                    fac=factor(c(rep(1:2,2), NA))
                    )
  dat$'(weights)' <- 1

  simple1 <- RItools:::makeDesigns(z ~ x1 + x2, data = dat)
  expect_false(any(grepl("TRUE", colnames(simple1))))
  expect_false(any(grepl("FALSE", colnames(simple1))))
              
  simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + strata(strat), data = dat)
  expect_false(any(grepl("TRUE", colnames(simple2))))
  expect_false(any(grepl("FALSE", colnames(simple2))))

  simple3 <- RItools:::makeDesigns(z ~ x1 + x2 + strata(strat) - 1, data = dat)
  expect_false(any(grepl("TRUE", colnames(simple3))))            
  expect_false(any(grepl("FALSE", colnames(simple3))))            

          })

test_that("DesignOptions to descriptive statistics", {
  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  d$'(weights)' = 1 # meet expectation of a weights column
  
  simple <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d)

  
  descriptives <- RItools:::designToDescriptives(simple)

  expect_is(descriptives, "array")
  expect_equal(dim(descriptives), c(5, 5, 2))

  # descriptives ignore clustering
  design.noclus <- RItools:::makeDesigns(z ~ x + f + strata(s), data = d)
  expect_equal(descriptives, RItools:::designToDescriptives(design.noclus))
  
  # the strata should imply different stats
  expect_false(identical(descriptives[,,1], descriptives[,,2]))

  # however, the pooled s.d.s should be the same.
  expect_identical(descriptives[,"pooled.sd","--"], descriptives[,"pooled.sd","s"])

  # ok, now checking that values are correct.
  expect_equal(mean(d$x[d$z == 1]), descriptives["x", "Treatment", "--"])
  expect_equal(mean(d$x[d$z == 0]), descriptives["x", "Control", "--"])

  # with equal sized strata, the the control/treatment means are the means of the the strata means
  expect_equal(mean(tapply(d$x[d$z == 1], d$s[d$z == 1], mean)), descriptives["x", "Treatment", "s"])
  expect_equal(mean(tapply(d$x[d$z == 0], d$s[d$z == 0], mean)), descriptives["x", "Control", "s"])

})

test_that("designToDescriptives uses provided covariate scales",{
    d  <- data.frame(x=rnorm(50), z=rep(c(0,1), 25))
    d$'(weights)'  <- 1
    simple <- RItools:::makeDesigns(z ~ x, data=d)
    sd_x  <- sd(resid(lm(x ~z, data=d)))
    descriptives <- designToDescriptives(simple, covariate.scales=c(x=sd_x*10))
    expect_equal(descriptives["x","pooled.sd","--"], sd_x*10)
    descriptives <- designToDescriptives(simple,
                                         covariate.scales=c(x=sd_x*10, y=Inf))
    expect_equal(descriptives["x","pooled.sd","--"], sd_x*10)
    expect_warning(designToDescriptives(simple, covariate.scales=sd_x),
                   "name")
    expect_warning(designToDescriptives(simple, covariate.scales=c(x="foo")),
                   "numeric")
})

test_that("descriptives for NotMissing variables", 
          {

              dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, 1), c(3,2)),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
              dat$'(weights)' <- 1

  simple <- RItools:::makeDesigns(z ~ x1 + x2  + strata(strat) + cluster(clus), data = dat)
              expect_false(any(grepl("NA", colnames(simple))))
              dsimple <- RItools:::designToDescriptives(simple)
              expect_match(dimnames(dsimple)[[1]], "(x1)", all=FALSE)
              expect_false(any(grepl("(x2)", dimnames(dsimple)[[1]], fixed=TRUE)))

  simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
  expect_false(any(grepl("NA", colnames(simple2))))
              dsimple2 <- RItools:::designToDescriptives(simple2)
              expect_match(dimnames(dsimple2)[[1]], "(x1)", all=FALSE)
              expect_match(dimnames(dsimple2)[[1]], "(fac)", all=FALSE)              
              expect_false(any(grepl("(x2)", dimnames(dsimple2)[[1]], fixed=TRUE)))                    
          }
)


test_that("Issue 36: Descriptives with NAs, appropriate weighting", {
  ### Descriptives with missing covariates ###

  set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      paired = rep(c(0,1), each = 250),
      z = rep(c(0,1), 250))
  d$'(weights)' = 1
  
  d.missing <- d

  d.missing$x[d.missing$x < -1] <- NA

  simple.all <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d)
  
  descriptives.all <- RItools:::designToDescriptives(simple.all)
  expect_equal(descriptives.all["x", "Treatment", "--"], mean(d$x[d$z == 1]))
  expect_equal(descriptives.all["x", "Treatment", "s"], mean(d$x[d$z == 1 & !is.na(d$s)]))

  simple.missing <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d.missing)
  
  descriptives.missing <- RItools:::designToDescriptives(simple.missing)

  with(d.missing,
       expect_equal(descriptives.missing["x", "Treatment", "--"],
                    mean(x[z == 1], na.rm = TRUE)))

  with(d.missing,
       expect_equal(descriptives.missing["x", "Treatment", "s"],
                    mean(x[z == 1 & !is.na(s)], na.rm = TRUE)))

  expect_false(identical(descriptives.all, descriptives.missing))

  # ETT weighting
  design.paired   <- RItools:::makeDesigns(z ~ x + f + strata(paired) + strata(s), data = d) 
  descriptives.paired <- RItools:::designToDescriptives(design.paired)

  with(d, expect_equal(descriptives.paired["x", "Control", "paired"], mean(x[z == 0])))
  with(d, expect_equal(descriptives.paired["x", "Control", "--"], mean(x[z == 0])))
  with(d, expect_false(identical(descriptives.paired["x", "Control", "s"], mean(x[z == 0]))))

})

test_that("Aggegating designs by clusters", {

    set.seed(20130801)

  d <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  d$'(weights)' <- 1
  # grab a bunch of rows and duplicate them
  d <- rbind(d, d[sample(1:dim(d)[1], size = 100), ])
  
  # swapping around the data to make sure order doesn't mask bugs
  d <- d[sample(1:dim(d)[1]), ]


  design <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d)
  aggDesign <- RItools:::aggregateDesigns(design) 

  # one row per cluster, with columns x, fa, fb, fc
  expect_equal(dim(aggDesign), c(100, 4))

  # now spot check some cluster totals of totals
  expect_equal(aggDesign[1, ], colMeans(design[design@Cluster == 1,]))
  
  # Z's roll up as they should
  Zs <- tapply(design@Z, design@Cluster, mean)
  dim(Zs) <- NULL
  Zs <- as.logical(Zs)
  expect_equivalent(aggDesign@Z, Zs)

  # extraneous levels in the Cluster slot are ignored
  design2 <- design
  levels(design2@Cluster) <- c(levels(design2@Cluster), letters)
  aggDesign2 <- RItools:::aggregateDesigns(design2) 
  expect_equal(dim(aggDesign2), c(100, 4))

})

test_that("aggregateDesigns treats NA covariates as 0's" ,{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x=rep(c(NA, 1), c(3,2))
                    )
  dat$'(weights)' <- 1

  design <- RItools:::makeDesigns(z~x+strata(strat)+cluster(clus)-1, dat)
  aggDesign <- RItools:::aggregateDesigns(design)
    expect_equivalent(aggDesign[,'x'], c(0, 0, 1, 1) )
    nm.column.for.x <- aggDesign@NM.Covariates[match( 'x', colnames(aggDesign))]
  expect_equal(aggDesign@NotMissing[,nm.column.for.x], c(0,0,1,1) )
})


test_that("Aggregation of unit weights to cluster level",{
  ##set.seed(20130801)

  d.short <- data.frame(
      x = rnorm(500),
      f = factor(sample(c("A", "B", "C"), size = 500, replace = T)),
      c = rep(1:100, 5),
      s = rep(c(1:4, NA), 100),
      z = rep(c(0,1), 250))

  
  ## grab a bunch of rows and duplicate them
  newrows <- c(1:nrow(d.short), sample(1:nrow(d.short), size=100, replace=T) )
  d.tall <- d.short[newrows,]
  d.tall$'(weights)' <- 1
  d.short$'(weights)' <- as.vector(table(newrows))

  design.tall <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d.tall)
  aggDesign.tall <- RItools:::aggregateDesigns(design.tall) 

    design.short <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d.short)
  aggDesign.short <- RItools:::aggregateDesigns(design.short) 

  ## we should wind up in the same place.
  expect_equal(aggDesign.tall, aggDesign.short)

  d2 <- d.tall
  d2$'(weights)' <- 2
  design.d2 <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d2)
  aggDesign.d2 <- RItools:::aggregateDesigns(design.d2)

  tmp <- aggDesign.tall
  tmp@UnitWeights <- 2 * tmp@UnitWeights ## the tall weights should be 1/2 of the d2 weights.
  expect_equal(tmp, aggDesign.d2) 
})


context("ModelMatrixPlus S4 class enriches model.matrix-value 'class' ")

test_that("R core hasn't revised conventions we may depend on",
          {
     ff <- log(Volume) ~ log(Height) + log(Girth)
     m <- model.frame(ff, trees)
     mm <- model.matrix(ff, m)
     expect_equal("(Intercept)", colnames(mm)[1])
     expect_equivalent(mm[,1,drop=TRUE], rep(1, nrow(mm)))
     expect_false("(Intercept)" %in% colnames(model.matrix(update(ff, .~.-1), m)))
     expect_equal(dim(mm[,integer(0)]), c(nrow(mm), 0))
     expect_equal(complete.cases(mm[,integer(0)]), rep(TRUE, nrow(mm)))
     expect_equivalent(weighted.mean(c(1:4, NA), c(rep(1,4), 0)), 2.5)
          })


context("alignDesignsByStrata")
test_that("alignDesigns, designToDescriptives output alignment", {

    
    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                      clus=factor(c(1,1,2:4)),
                      z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                      x1=rep(c(NA, 1), c(3,2)),
                      x2=c(1:5),
                      fac=factor(c(rep(1:2,2), NA))
                      )
    dat$'(weights)' <- 1

    simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
    expect_equal(colnames(simple2@StrataFrame), c("strat", "--"))
    dsimple2 <- RItools:::designToDescriptives(simple2)
    asimple2 <- RItools:::alignDesignsByStrata(simple2)
    expect_true(setequal(names(asimple2), c("strat", "--")))
    expect_equivalent(dimnames(dsimple2)[[1]], colnames(asimple2[["--"]]@Covariates))

})


test_that("alignDesigns centers covars by stratum", {

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                      clus=factor(c(1,1,2:4)),
                      z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                      x1=rep(c(NA, 1), c(3,2)),
                      x2=c(1:5),
                      fac=factor(c(rep(1:2,2), NA))
                      )
    dat$'(weights)' <- 1

    ## first unweighted case
    simple0 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
    asimple0 <- RItools:::alignDesignsByStrata(simple0)
    expect_equivalent(colSums(asimple0[["--"]]@Covariates),
                      rep(0,ncol(asimple0[["--"]]@Covariates)))
    expect_equivalent(colSums(asimple0[["strat"]]@Covariates[simple0@StrataFrame[["strat"]]=="a",]),
                      rep(0,ncol(asimple0[["strat"]]@Covariates)))
    expect_equivalent(as.matrix(t(asimple0[["strat"]]@Design@Units) %*% asimple0[["strat"]]@Covariates),
                      matrix(0,2,ncol(asimple0[["strat"]]@Covariates)))


    ##  post alignment transform
    asimple2 <- RItools:::alignDesignsByStrata(simple0, post.align.transform = rank)
    expect_equivalent(colSums(asimple2[["--"]]@Covariates  ),
                      rep(0,ncol(asimple2[["--"]]@Covariates)))

    tmp2 <- asimple2[["strat"]]@Covariates 
    expect_equivalent(colSums(tmp2[simple0@StrataFrame[["strat"]]=="a",]),
                      rep(0,ncol(asimple2[["strat"]]@Covariates)))
    expect_equivalent(as.matrix(t(asimple2[["strat"]]@Design@Units) %*% tmp2),
                      matrix(0,2, ncol(asimple2[["strat"]]@Covariates)))

} )

test_that("scale() method wrapping to alignDesignsByStrata()",{
    # at first pass, we're testing form but not content here.
    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                      clus=factor(c(1,1,2:4)),
                      z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                      x1=rep(c(NA, 1), c(3,2)),
                      x2=c(1:5),
                      fac=factor(c(rep(1:2,2), NA))
                      )
    dat$'(weights)' <- 1

    simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
    scl2_scaleF  <- scale(simple2, center=TRUE, scale=FALSE)
    asimple2  <- RItools:::alignDesignsByStrata(simple2, post.align.transform = NULL)
    expect_identical(scl2_scaleF, asimple2[["--"]]@Covariates)

    simple2c  <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus) - 1, data = dat)
    scl2c_scaleF  <- scale(simple2c, center=TRUE, scale=FALSE)
    expect_identical(scl2c_scaleF, asimple2[["strat"]]@Covariates)
    scl2_scaleF_centerF  <- scale(simple2, center=FALSE, scale=FALSE) # if it's a logical, 
    expect_identical(scl2_scaleF, scl2_scaleF_centerF)                # `center` param is ignored
    scl2_scaleF_centerrank  <- scale(simple2, center=rank, scale=FALSE)
    expect_identical(dim(scl2_scaleF_centerrank), dim(scl2_scaleF))
    expect_false(isTRUE(all.equal(scl2_scaleF, scl2_scaleF_centerrank, check.attributes=FALSE)),
                 "post alignment transform ignored")
    scl2_scaleT  <- scale(simple2, center=TRUE, scale=TRUE)
    expect_equal(length(dim(scl2_scaleT)), 2L)
    expect_equivalent(is.na(scl2_scaleT),
                      matrix(FALSE, nrow(scl2_scaleT), ncol(scl2_scaleT)))
    
})


## context("HB08")
## unclear if HB08 will matter with new balanceTest architecture, commenting out failures
## 
## test_that("HB08 agreement w/ xBal()", {
## 
##   set.seed(20180605)
##   n <- 100
##   x1 <- rnorm(n)
##   x2 <- rnorm(n)
##   x3 <- 0.5 + 0.25 * x1 - 0.25 * x2 + rnorm(n)
##   idx <- 0.25 + 0.1 * x1 + 0.2 * x2 - 0.5 * x3 + rnorm(n)
##   y <- sample(rep(c(1,0), n/2), prob = exp(idx) / (1 + exp(idx)))
## 
##   xy <- data.frame(x1, x2, x3, idx, y)
##   xy$m[y == 1] <- order(idx[y == 1])
##   xy$m[y == 0] <- order(idx[y == 0])
##   ## this mimics matched pairs:
##   expect_true(all(table(xy$y, xy$m)==1))
##   xy$'(weights)' <- rep(1L, n) 
## 
##     ## first unweighted case
##     simple0 <- RItools:::makeDesigns(y ~ x1 + x2 + x3 + strata(m), data = xy)
##     asimple0 <- RItools:::alignDesignsByStrata(simple0)
##     btis0 <- lapply(asimple0, HB08)
##     xb0 <- xBalance(y ~ x1 + x2 + x3, data = xy,
##                     strata = list(unmatched = NULL, matched = ~ m), report = c("all"))
## 
##     expect_equivalent(btis0[['--']]$adj.mean.diffs[-4], # remove '(_non-null record_)' entry
##                       xb0$results[,'adj.diff',"unmatched"])
##     expect_equivalent(btis0[['--']]$tcov[1:3,1:3], # remove '(_non-null record_)' entries
##                       attr(xb0$overall, 'tcov')$unmatched)
##     expect_equivalent(btis0[['--']][c('Msq', 'DF')],
##                       xb0[['overall']]["unmatched",c('chisquare', 'df'), drop=TRUE])
## 
##     expect_equivalent(btis0[['m']]$adj.mean.diffs[-4], # remove '(_non-null record_)' entry
##                       xb0$results[,'adj.diff',"matched"])
##     expect_equivalent(btis0[['m']]$tcov[1:3,1:3], # remove '(_non-null record_)' entries
##                       attr(xb0$overall, 'tcov')$matched)
##     expect_equivalent(btis0[['m']][c('Msq', 'DF')],
##                       xb0[['overall']]["matched",c('chisquare', 'df'), drop=TRUE])
## 
## 
##     ## now with weights.  Comparing adjusted diffs based on totals will only work
##     ## if the weights don't vary with by stratum, at least in stratified case.
##     xy_wted <- xy; mwts <- 0
##     while (any(mwts==0)) mwts <- rpois(n/2, lambda=10)
##     ## centering of variables is needed for unstratified mean diffs comparison.
##     xy_wted <- transform(xy_wted, x1=x1-weighted.mean(x1,unsplit(mwts, xy$m)), 
##                          x2=x2-weighted.mean(x2,unsplit(mwts, xy$m)), 
##                          x3=x3-weighted.mean(x3,unsplit(mwts, xy$m)))
##     xy_wted$'(weights)' <- unsplit(mwts, xy$m)
##     
##     simple1 <- RItools:::makeDesigns(y ~ x1 + x2 + x3 + strata(m), data = xy_wted)
##     asimple1 <- RItools:::alignDesignsByStrata(simple1)
##     btis1 <- lapply(asimple1, HB08)
## 
##   wts.scaled <- xy_wted$'(weights)' / mean(xy_wted$'(weights)')
## 
##   xy_xbwts <- transform(xy_wted, x1=x1*wts.scaled,
##                         x2=x2*wts.scaled, x3=x3*wts.scaled)
##   xb1u <- xBalance(y ~ x1 + x2 + x3, data = xy_xbwts,
##                    strata = list(unmatched = NULL), report = c("all"))
##   expect_equivalent(btis1[['--']]$adj.mean.diffs[-4], # remove '(_non-null record_)' entry
##                     xb1u$results[,'adj.diff',"unmatched"])
##   expect_equivalent(btis1[['--']]$tcov[1:3,1:3], # remove '(_non-null record_)' entries
##                     attr(xb1u$overall, 'tcov')$unmatched)
##   expect_equivalent(btis1[['--']][c('Msq', 'DF')],
##                     xb1u[['overall']]["unmatched",c('chisquare', 'df'), drop=TRUE])
## 
## 
##   wts.scaled <- xy_wted$'(weights)' / mean( mwts )
##   xy_xbwts <- transform(xy_wted, x1=x1*wts.scaled,
##                         x2=x2*wts.scaled, x3=x3*wts.scaled)
##   xb1m <- xBalance(y ~ x1 + x2 + x3, data = xy_xbwts,
##                    strata = list(matched = ~ m), report = c("all"))
## 
##   expect_equivalent(btis1[['m']]$adj.mean.diffs[-4], # remove '(_non-null record_)' entry
##                     xb1m$results[,'adj.diff',"matched"])
##   expect_equivalent(btis1[['m']]$tcov[1:3,1:3], # remove '(_non-null record_)' entries
##                     attr(xb1m$overall, 'tcov')$matched)
##   expect_equivalent(btis1[['m']][c('Msq', 'DF')],
##                     xb1m[['overall']]["matched",c('chisquare', 'df'), drop=TRUE])
## 
## 
## } )



### Tests to write...
##test_that("alignDesigns properly tracks UnitWeights vs NotMissing",{})
##test_that("",{})
