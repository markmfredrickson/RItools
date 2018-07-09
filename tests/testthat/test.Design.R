################################################################################
# Tests for DesignOptions objects
################################################################################

library("testthat")

context("DesignMatrix S4 class carries per-covariate missingness info")

test_that("If no missing data, then NotMissing is a matrix w n rows and 0 cols",{

    d <- data.frame(
      id          = 1:500,
      x           = rnorm(500),
      cluster     = rep(1:100, 5),
      strata.good = rep(c(1:4, NA), 100),
      strata.bad  = sample(1:100, size = 500, replace = T),
      z.good      = rep(c(0,1), 250),
      z.bad       = sample(c(0,1), size = 500, replace = T))
  d$'(weights)' = 1 # meet expectation of a weights column
  simple <- RItools:::design_matrix(z.good ~ x, data = d)
  expect_equivalent(dim(simple@NotMissing), c(500,1))

})
test_that("Missingness gets passed through in Covariates, recorded in NotMissing",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, 1), c(3,2)),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
              dat$'(weights)' <- 1
              datmf <- model.frame(z ~ x1 + x2 + fac, dat, na.action = na.pass) 
              simple2 <- RItools:::design_matrix(z ~ x1 + x2 + fac, data = datmf)
              expect_equivalent(ncol(simple2@NotMissing), 3)
              expect_equivalent(colnames(simple2@NotMissing), c("_any Xs recorded_", "x1", "fac"))
              
})
test_that("lookup tables OK, even w/ complex & multi-column terms",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, 1), c(3,2)),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
    dat$'(weights)' <- 1
    datmf <- model.frame(z ~ x1 + x2 + fac, dat, na.action = na.pass)
    simple2 <- RItools:::design_matrix(z ~ x1 + x2 + fac, data=datmf)
    expect_equal(simple2@OriginalVariables, 1:3)
    expect_equal(simple2@TermLabels, c( "x1", "x2", "fac"))
    expect_equal(simple2@NM.Covariates, c(2,0,3))
    expect_equal(simple2@NM.terms, c(2,0,3))

    ## check that complex term don't spell trouble in themselves
    datmf <- model.frame(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = dat,
                         na.action = na.pass)
    simple3 <- RItools:::design_matrix(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = datmf)
    expect_equal(simple3@OriginalVariables, 1:3)
    expect_equal(simple3@TermLabels, c("x1", "cut(x2, c(0, 3, 6))", "fac"))
    expect_equal(simple3@NM.Covariates, c(2,0,3))
    expect_equal(simple3@NM.terms, c(2,0,3))

    ## now try a complex term that actually expands to multiple columns
    datmf <- model.frame(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = dat,
                         na.action = na.pass)
    simple4 <- RItools:::design_matrix(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = datmf,
                                       contrasts=list("cut(x2, c(0, 3, 6))"=diag(2)))
    expect_equal(simple4@OriginalVariables, c(1,2,2,3))
    expect_equal(simple4@TermLabels, c("x1", "cut(x2, c(0, 3, 6))", "fac"))
    expect_equal(simple4@NM.Covariates, c(2,0,0,3))
    expect_equal(simple4@NM.terms, c(2,0,3))

    ## now put NAs in the multi-column complex term
    datmf <- model.frame(z ~ x2 + cut(x1, c(0,3,6)) + fac, data = dat,
                         na.action = na.pass)
    simple5 <- RItools:::design_matrix(z ~ x2 + cut(x1, c(0,3,6)) + fac, data = datmf,
                                       contrasts=list("cut(x1, c(0, 3, 6))"=diag(2)))
    expect_equal(simple5@OriginalVariables, c(1,2,2,3))
    expect_equal(simple5@TermLabels, c("x2", "cut(x1, c(0, 3, 6))", "fac"))
    expect_equal(simple5@NM.Covariates, c(0,2,2,3))
    expect_equal(simple5@NM.terms, c(0,2,3))
 
})

test_that("Issue #76: Using I() in formulas", {
  x <- data.frame(x = rnorm(10), y = rnorm(10), z = rbinom(10, size = 1, p = 1/3))
  x$"(weights)" <- 1
  d <- makeDesigns(z ~ I(x * sin(y)), data = x)
  expect_s4_class(d, "DesignOptions")
})

test_that("Issue #86: makeDesigns finds variables outside data arg",{
    data(nuclearplants)
    foo <- nuclearplants$pt
    nuclearplants$"(weights)" <- 1
    d <- makeDesigns(pr ~ cost + strata(foo), data=nuclearplants)
      expect_s4_class(d, "DesignOptions")
})

test_that("Duplicated missingness patterns handled appropriately",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=c(rep(NA, 3), 1:2),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
    dat$'(weights)' <- 1
    datmf <- model.frame(z ~ x1 + I(x1^2) + fac, data=dat, na.action = na.pass)
    simple2 <- RItools:::design_matrix(z ~ x1 + I(x1^2) + fac, data = datmf)
    expect_equal(simple2@OriginalVariables, 1:3)
    expect_equal(simple2@TermLabels, c("x1", "I(x1^2)", "fac"))
    expect_equal(simple2@NM.Covariates, 1+c(1,1,2))
    expect_equal(simple2@NM.terms, 1+c(1,1,2))

    ## If exactly two terms have missing data but in the same pattern, then
    ## NotMissing is a matrix w/ n rows and 1 col.
    datmf <- model.frame(z ~ x1 + I(x1^2), data=dat, na.action = na.pass)
    simple3 <- RItools:::design_matrix(z ~ x1 + I(x1^2), data = datmf)
    expect_equal(simple3@OriginalVariables, 1:2)
    expect_equal(simple3@TermLabels, c("x1", "I(x1^2)"))
    expect_equal(ncol(simple3@NotMissing), 1)
    expect_equal(simple3@NM.Covariates, c(1,1))
    expect_equal(simple3@NM.terms, c(1,1))

})

test_that("All-Xes missingness |-> NotMissing col '_any Xs recorded_'",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=c(rep(NA, 4), 2),
                    x2=c(1:3, NA, 5),
                    fac=factor(c(NA, 1:2, NA, 2))
                    )
    dat$'(weights)' <- 1
    datmf <- model.frame(z ~ x1 + I(x1^2) + fac, data=dat, na.action = na.pass)
    simple6 <- RItools:::design_matrix(z ~ x1 + I(x1^2) + fac, data = datmf)
    expect_equal(simple6@OriginalVariables, 1:3)
    expect_equal(simple6@TermLabels, c("x1", "I(x1^2)", "fac"))
    expect_equal(colnames(simple6@NotMissing)[1], "_any Xs recorded_")
    expect_equal(simple6@NM.Covariates, c(2,2,1))
    expect_equal(simple6@NM.terms, c(2,2,1))    

    simple7 <- makeDesigns(z ~ x1 + I(x1^2) + fac + 0 + strata(strat) + cluster(clus),
                           dat)
    simple7 <- as(simple7, "StratumWeightedDesignOptions")
    simple7@Sweights <-
        RItools:::DesignWeights(simple7, 
                                RItools:::effectOfTreatmentOnTreated)
    ## As writing of this test, DesignWeights() expects only pre-aggregated designs,
    ## and infers treatment:control ratios from the numbers of elements in in each
    ## of condition (by strata), not the number of clusters. Thus in this example
    ## it should believe that the ETT weights are proportional to 2 for stratum a,
    ## 1 for stratum b:
    expect_equivalent(simple7@Sweights$strat$sweights, (2:1)/3)
    ## key point: had missingness of each of unit 4's covariates tricked it into ignoring
    ## that observation, then the ETT weight for stratum b would have been proportional to 0,
    ## not 1.

} )

test_that("All-Xes missingness noted in descriptives", {
    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=c(rep(NA, 4), 2),
                    x2=c(1:3, NA, 5),
                    fac=factor(c(NA, 1:2, NA, 2))
                    )
    dat$'(weights)' <- 1

    simple5 <- makeDesigns(z ~ x1 + I(x1^2) + fac + cluster(clus), dat)
    descr <- RItools::designToDescriptives(simple5)
    expect_equivalent(descr['(_any Xs recorded_)' ,1:2,1], c(1, 1/3)) # not just 1s all around
})

context("DesignOptions objects")

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
  expect_equivalent(simple@Covariates[, "x"], d$x)
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
  expect_false(any(grepl("TRUE", colnames(simple@Covariates))))
  expect_false(any(grepl("FALSE", colnames(simple@Covariates))))

  simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
  expect_match(colnames(simple2@NotMissing), "x1", all=FALSE)
  expect_false(any(grepl("x2", colnames(simple2@NotMissing))))
  expect_match(colnames(simple2@NotMissing), "fac", all=FALSE)
  expect_false(any(grepl("TRUE", colnames(simple2@Covariates))))            
  expect_false(any(grepl("FALSE", colnames(simple2@Covariates))))

  simple3 <- RItools:::makeDesigns(z ~ x1 + x3 + fac+ strata(strat) + cluster(clus), data = dat)
  expect_match(colnames(simple3@NotMissing), "x1", all=FALSE)
  expect_match(colnames(simple3@NotMissing), "x3", all=FALSE)
  expect_match(colnames(simple3@NotMissing), "fac", all=FALSE)
  expect_false(any(grepl("TRUE", colnames(simple3@Covariates))))            
  expect_false(any(grepl("FALSE", colnames(simple3@Covariates))))
              
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
  expect_false(any(grepl("TRUE", colnames(simple1@Covariates))))
  expect_false(any(grepl("FALSE", colnames(simple1@Covariates))))
              
  simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + strata(strat), data = dat)
  expect_false(any(grepl("TRUE", colnames(simple2@Covariates))))
  expect_false(any(grepl("FALSE", colnames(simple2@Covariates))))

  simple3 <- RItools:::makeDesigns(z ~ x1 + x2 + strata(strat) - 1, data = dat)
  expect_false(any(grepl("TRUE", colnames(simple3@Covariates))))            
  expect_false(any(grepl("FALSE", colnames(simple3@Covariates))))            

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
  
  # ok, now checking that values are correct.
  expect_equal(mean(d$x[d$z == 1]), descriptives["x", "Treatment", "--"])
  expect_equal(mean(d$x[d$z == 0]), descriptives["x", "Control", "--"])

  # with equal sized strata, the the control/treatment means are the means of the the strata means
  expect_equal(mean(tapply(d$x[d$z == 1], d$s[d$z == 1], mean)), descriptives["x", "Treatment", "s"])
  expect_equal(mean(tapply(d$x[d$z == 0], d$s[d$z == 0], mean)), descriptives["x", "Control", "s"])

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
              expect_false(any(grepl("NA", colnames(simple@Covariates))))
              dsimple <- RItools:::designToDescriptives(simple)
              expect_match(dimnames(dsimple)[[1]], "(x1)", all=FALSE)
              expect_false(any(grepl("(x2)", dimnames(dsimple)[[1]], fixed=TRUE)))

  simple2 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat)
  expect_false(any(grepl("NA", colnames(simple2@Covariates))))
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
  expect_equal(dim(aggDesign@Covariates), c(100, 4))

  # now spot check some cluster totals of totals
  expect_equal(aggDesign@Covariates[1, ], colMeans(design@Covariates[design@Cluster == 1,]))
  
  # Z's roll up as they should
  Zs <- tapply(design@Z, design@Cluster, mean)
  dim(Zs) <- NULL
  Zs <- as.logical(Zs)
  expect_equivalent(aggDesign@Z, Zs)

  # extraneous levels in the Cluster slot are ignored
  design2 <- design
  levels(design2@Cluster) <- c(levels(design2@Cluster), letters)
  aggDesign2 <- RItools:::aggregateDesigns(design2) 
  expect_equal(dim(aggDesign2@Covariates), c(100, 4))

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
    expect_equivalent(aggDesign@Covariates[,'x'], c(0, 0, 1, 1) )
    nm.column.for.x <- aggDesign@NM.Covariates[match( 'x', colnames(aggDesign@Covariates))]
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

  expect_equal(2*aggDesign.tall@UnitWeights,aggDesign.d2@UnitWeights)
  expect_equal(aggDesign.tall@NotMissing, aggDesign.d2@NotMissing)
  expect_equal(aggDesign.tall@Covariates, aggDesign.d2@Covariates) 
          })


context("DesignMatrix S4 class enriches model.matrix-value 'class' ")

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

test_that("Model matrix material is properly formed",
          {
              
     ff <- log(Volume) ~ log(Height) + log(Girth)
     DM0 <- design_matrix(ff, trees, remove.intercept=FALSE)
     expect_is(DM0, "DesignMatrix")
     m <- model.frame(ff, trees)
     expect_equivalent(model.matrix(ff, m), as.matrix(DM0))
     fff <- update(ff, .~.-1)
     trees2 <- trees
     trees2[1, "Volume"] <- NA # LHS variable, shouldn't affect anything
     expect_equivalent(model.matrix(fff, m), as.matrix(design_matrix(fff, model.frame(fff, trees2, na.action = na.pass))))
     trees2[1, "Height"] <- NA # RHS variable, but still shouldn't cause rows to be dropped
     expect_equal(dim(model.matrix(fff, m)), dim(as.matrix(design_matrix(fff,
model.frame(fff, trees2, na.action = na.pass)
                                                                                ))))


     ## specified contrasts
     dd <- data.frame(a = gl(3,4), b = gl(4,1,12)) # balanced 2-way
     expect_equal(model.matrix(~ a + b, dd),
                  as.matrix(design_matrix(~ a + b, dd, remove.intercept=FALSE)))
     expect_equal(model.matrix(~ a + b-1, dd),
                  as.matrix(design_matrix(~ a + b-1, dd)))
     expect_equal(model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum")),
                  as.matrix(design_matrix(~ a + b, dd,
                                          contrasts = list(a = "contr.sum"), remove.intercept=FALSE))
                  )
     expect_equal(model.matrix(~ a + b, dd,
                               contrasts = list(a = "contr.sum", b = "contr.poly")),
                  as.matrix(design_matrix(~ a + b, dd,
                                          contrasts =
                                              list(a = "contr.sum", b = "contr.poly"),
                                          remove.intercept=FALSE))
                  )
          }
          )

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
    simple2 <-   as(simple2, "StratumWeightedDesignOptions")
    simple2@Sweights <- RItools:::DesignWeights(simple2, # Have to aggregate 1st to figure stratum weights
                                                RItools:::effectOfTreatmentOnTreated)
    expect_true(setequal(names(simple2@Sweights), c("strat", "--")))
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
    simple0 <-   as(simple0, "StratumWeightedDesignOptions")
    simple0@Sweights <- RItools:::DesignWeights(simple0, # Placeholder strat weights, shouldn't affect 
                                                RItools:::effectOfTreatmentOnTreated) # this test
    asimple0 <- RItools:::alignDesignsByStrata(simple0)
    expect_equivalent(colSums(asimple0[["--"]]@Covariates),
                      rep(0,ncol(asimple0[["--"]]@Covariates)))
    expect_equivalent(colSums(asimple0[["strat"]]@Covariates[asimple0[["strat"]]@StrataFactor=="a",]),
                      rep(0,ncol(asimple0[["strat"]]@Covariates)))
    expect_equivalent(colSums(asimple0[["strat"]]@Covariates[asimple0[["strat"]]@StrataFactor=="b",]),
                      rep(0,ncol(asimple0[["strat"]]@Covariates)))

    ## now with weights
    dat1 <- dat
    dat1$'(weights)' <- rpois(nrow(dat1), lambda=10)
    while (any(dat1$'(weights)'==0)) dat1$'(weights)' <- rpois(nrow(dat1), lambda=10)
    
    simple1 <- RItools:::makeDesigns(z ~ x1 + x2 + fac+ strata(strat) + cluster(clus), data = dat1)
    simple1 <-   as(simple1, "StratumWeightedDesignOptions")
    simple1@Sweights <- RItools:::DesignWeights(simple1, # Placeholder strat weights, shouldn't affect 
                                                RItools:::effectOfTreatmentOnTreated) # this test
    asimple1 <- RItools:::alignDesignsByStrata(simple1)
    expect_equivalent(colSums(asimple1[["--"]]@Covariates *
                              asimple1[["--"]]@UnitWeights ),
                      rep(0,ncol(asimple1[["--"]]@Covariates)))

    tmp1 <- asimple1[["strat"]]@Covariates * asimple1[["strat"]]@UnitWeights 
    expect_equivalent(colSums(tmp1[asimple1[["strat"]]@StrataFactor=="a",]),
                      rep(0,ncol(asimple1[["strat"]]@Covariates)))
    expect_equivalent(colSums(tmp1[asimple1[["strat"]]@StrataFactor=="b",]),
                      rep(0,ncol(asimple1[["strat"]]@Covariates)))

    ## now with weights, post alignment transform
    asimple2 <- RItools:::alignDesignsByStrata(simple1, post.align.transform = rank)
    expect_equivalent(colSums(asimple2[["--"]]@Covariates *
                              asimple2[["--"]]@UnitWeights ),
                      rep(0,ncol(asimple2[["--"]]@Covariates)))

    tmp2 <- asimple2[["strat"]]@Covariates * asimple2[["strat"]]@UnitWeights 
    expect_equivalent(colSums(tmp2[asimple2[["strat"]]@StrataFactor=="a",]),
                      rep(0,ncol(asimple2[["strat"]]@Covariates)))
    expect_equivalent(colSums(tmp2[asimple2[["strat"]]@StrataFactor=="b",]),
                      rep(0,ncol(asimple2[["strat"]]@Covariates)))

} )


test_that("Issue #89: Proper strata weights", {

  set.seed(20180208)

  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- 0.5 + 0.25 * x1 - 0.25 * x2 + rnorm(n)
  idx <- 0.25 + 0.1 * x1 + 0.2 * x2 - 0.5 * x3 + rnorm(n)
  y <- sample(rep(c(1,0), n/2), prob = exp(idx) / (1 + exp(idx)))

  xy <- data.frame(x1, x2, x3, idx, y)
  xy$m[y == 1] <- order(idx[y == 1])
  xy$m[y == 0] <- order(idx[y == 0])
  xy$"(weights)" <- 1

  xy.wts <- xy
  xy.wts$"(weights)" <- (1 + exp(idx)) / exp(idx) # inverse propensity score weights

  design.nowts <- RItools:::makeDesigns(y ~ x1 + x2 + x3 + strata(m), data  = xy)
  design.wts <- RItools:::makeDesigns(y ~ x1 + x2 + x3 + strata(m), data = xy.wts)

  ## ETT weights are determined by assignment probabilities, not element counts 
  ## or cluster masses. Accordingly presence/absence of unit weights shouldn't matter
  ## for their sweights.  (But since they do affect h_b * m-bar_b, the corresponding
  ## wtratio's will be affected.)
  ett.nowts <- RItools:::DesignWeights(design.nowts, stratum.weights=effectOfTreatmentOnTreated)
  ett.wts <- RItools:::DesignWeights(design.wts, stratum.weights=effectOfTreatmentOnTreated)

  expect_equal(ett.wts$m$sweights, ett.nowts$m$sweights)

  ## With a single stratum, sweights has to be 1, since it's normalized.
  ## wtratio is its ratio with h_b * m-bar_b, thus will generally be much 
  ## less than 1.  Check this:
  h <- with(xy, 1/(1/sum(y) + 1/sum(!y)))
  expect_equal(ett.nowts[['--']][,'wtratio'], 1/h)
  h <- with(xy.wts, 1/(1/sum(y) + 1/sum(!y)))
  expect_equal(ett.wts[['--']][,'wtratio'], 1/(h*mean(xy.wts$"(weights)")))
  
  
  ## split up into strata, use harmonic strata weights
  ## again unit weights shouldn't enter into this, although they
  ## would affect harmonic_times_mean_weight
  dw.nowts <- RItools:::DesignWeights(design.nowts, stratum.weights=harmonic)
  dw.wts <- RItools:::DesignWeights(design.wts, stratum.weights=harmonic)

  expect_equal(dw.wts$m$sweights, dw.nowts$m$sweights)

  ## in this example by-stratum harmonic mean cluster counts are always 1 --
  
  expect_equivalent(as.vector(dw.wts$m$sweights),
               rep(1, nlevels(survival::strata(xy.wts$m)))/
                 nlevels(survival::strata(xy.wts$m))
               )
  ## -- so we can check the calculation of the mean cluster mass factor as
  ## follows. 
  dw.wts2 <- RItools:::DesignWeights(design.wts)
  clus_mean_weights <- tapply(xy.wts$"(weights)", xy.wts$m, mean)
  expect_equivalent(dw.wts2$m$sweights, clus_mean_weights/sum(clus_mean_weights))

})

context("alignedToInferentials")


test_that("alignedToInferentials agreement w/ xBal()", {

  set.seed(20180605)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- 0.5 + 0.25 * x1 - 0.25 * x2 + rnorm(n)
  idx <- 0.25 + 0.1 * x1 + 0.2 * x2 - 0.5 * x3 + rnorm(n)
  y <- sample(rep(c(1,0), n/2), prob = exp(idx) / (1 + exp(idx)))

  xy <- data.frame(x1, x2, x3, idx, y)
  xy$m[y == 1] <- order(idx[y == 1])
  xy$m[y == 0] <- order(idx[y == 0])
  ## this mimics matched pairs:
  expect_true(all(table(xy$y, xy$m)==1))
  xy$'(weights)' <- rep(1L, n) 

    ## first unweighted case
    simple0 <- RItools:::makeDesigns(y ~ x1 + x2 + x3 + strata(m), data = xy)
    simple0 <-   as(simple0, "StratumWeightedDesignOptions")
    simple0@Sweights <- RItools:::DesignWeights(simple0) # this test
    asimple0 <- RItools:::alignDesignsByStrata(simple0)
    btis0 <- lapply(asimple0, alignedToInferentials)
    xb0 <- xBalance(y ~ x1 + x2 + x3, data = xy,
                    strata = list(unmatched = NULL, matched = ~ m), report = c("all"))

    expect_equivalent(btis0[['--']]$adj.mean.diffs[-4], # remove '(_any Xs recorded_)' entry
                      xb0$results[,'adj.diff',"unmatched"])
    expect_equivalent(btis0[['--']]$tcov[1:3,1:3], # remove '(_any Xs recorded_)' entries
                      attr(xb0$overall, 'tcov')$unmatched)
    expect_equivalent(btis0[['--']][c('csq', 'DF')],
                      xb0[['overall']]["unmatched",c('chisquare', 'df'), drop=TRUE])

    expect_equivalent(btis0[['m']]$adj.mean.diffs[-4], # remove '(_any Xs recorded_)' entry
                      xb0$results[,'adj.diff',"matched"])
    expect_equivalent(btis0[['m']]$tcov[1:3,1:3], # remove '(_any Xs recorded_)' entries
                      attr(xb0$overall, 'tcov')$matched)
    expect_equivalent(btis0[['m']][c('csq', 'DF')],
                      xb0[['overall']]["matched",c('chisquare', 'df'), drop=TRUE])

    ## now with weights.  Comparing adjusted diffs based on totals will only work
    ## if the weights don't vary with by stratum, at least in stratified case.
    xy_wted <- xy; mwts <- 0
    while (any(mwts==0)) mwts <- rpois(n/2, lambda=10)
    ## centering of variables is needed for unstratified mean diffs comparison.
    xy_wted <- transform(xy_wted, x1=x1-weighted.mean(x1,unsplit(mwts, xy$m)), 
                         x2=x2-weighted.mean(x2,unsplit(mwts, xy$m)), 
                         x3=x3-weighted.mean(x3,unsplit(mwts, xy$m)))
    xy_wted$'(weights)' <- unsplit(mwts, xy$m)
    
    simple1 <- RItools:::makeDesigns(y ~ x1 + x2 + x3 + strata(m), data = xy_wted)
    simple1 <-   as(simple1, "StratumWeightedDesignOptions")
    simple1@Sweights <- RItools:::DesignWeights(simple1) # this test
    asimple1 <- RItools:::alignDesignsByStrata(simple1)
    btis1 <- lapply(asimple1, alignedToInferentials)

  wts.scaled <- xy_wted$'(weights)' / mean(xy_wted$'(weights)')

  xy_xbwts <- transform(xy_wted, x1=x1*wts.scaled,
                        x2=x2*wts.scaled, x3=x3*wts.scaled)
  xb1u <- xBalance(y ~ x1 + x2 + x3, data = xy_xbwts,
                   strata = list(unmatched = NULL), report = c("all"))
  expect_equivalent(btis1[['--']]$adj.mean.diffs[-4], # remove '(_any Xs recorded_)' entry
                    xb1u$results[,'adj.diff',"unmatched"])
  expect_equivalent(btis1[['--']]$tcov[1:3,1:3], # remove '(_any Xs recorded_)' entries
                    attr(xb1u$overall, 'tcov')$unmatched)
  expect_equivalent(btis1[['--']][c('csq', 'DF')],
                    xb1u[['overall']]["unmatched",c('chisquare', 'df'), drop=TRUE])


  wts.scaled <- xy_wted$'(weights)' / mean( mwts )
  xy_xbwts <- transform(xy_wted, x1=x1*wts.scaled,
                        x2=x2*wts.scaled, x3=x3*wts.scaled)
  xb1m <- xBalance(y ~ x1 + x2 + x3, data = xy_xbwts,
                   strata = list(matched = ~ m), report = c("all"))

  expect_equivalent(btis1[['m']]$adj.mean.diffs[-4], # remove '(_any Xs recorded_)' entry
                    xb1m$results[,'adj.diff',"matched"])
  expect_equivalent(btis1[['m']]$tcov[1:3,1:3], # remove '(_any Xs recorded_)' entries
                    attr(xb1m$overall, 'tcov')$matched)
  expect_equivalent(btis1[['m']][c('csq', 'DF')],
                    xb1m[['overall']]["matched",c('chisquare', 'df'), drop=TRUE])


} )

### Tests to write...
##test_that("alignDesigns properly tracks UnitWeights vs NotMissing",{})
##test_that("",{})
