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
test_that("Missingingess gets passed through in Covariates, recorded in NotMissing",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, 1), c(3,2)),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
              dat$'(weights)' <- 1
              simple2 <- RItools:::design_matrix(z ~ x1 + x2 + fac, data=dat)
              expect_equivalent(ncol(simple2@NotMissing), 3)
              expect_equivalent(colnames(simple2@NotMissing), c("element weight", "x1", "fac"))
              
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
    simple2 <- RItools:::design_matrix(z ~ x1 + x2 + fac, data=dat)
    expect_equal(simple2@OriginalVariables, 1:3)
    expect_equal(simple2@TermLabels, c("x1", "x2", "fac"))
    expect_equal(simple2@NM.Covariates, 1+c(1,0,2))
    expect_equal(simple2@NM.terms, 1+c(1,0,2))
    ## check that complex term don't spell trouble in themselves
    simple3 <- RItools:::design_matrix(z ~ x1 + cut(x2, c(0,3,6)) + fac, data=dat)
    expect_equal(simple3@OriginalVariables, 1:3)
    expect_equal(simple3@TermLabels, c("x1", "cut(x2, c(0, 3, 6))", "fac"))
    expect_equal(simple3@NM.Covariates, 1+c(1,0,2))
    expect_equal(simple3@NM.terms, 1+c(1,0,2))
    ## now try a complex term that actually expands to multiple columns
    simple4 <- RItools:::design_matrix(z ~ x1 + cut(x2, c(0,3,6)) + fac, data=dat,
                                       contrasts=list("cut(x2, c(0, 3, 6))"=diag(2)))
    expect_equal(simple4@OriginalVariables, c(1,2,2,3))
    expect_equal(simple4@TermLabels, c("x1", "cut(x2, c(0, 3, 6))", "fac"))
    expect_equal(simple4@NM.Covariates, 1+c(1,0,0,2))
    expect_equal(simple4@NM.terms, 1+c(1,0,2))
    ## now put NAs in the multi-column complex term
    simple5 <- RItools:::design_matrix(z ~ x2 + cut(x1, c(0,3,6)) + fac, data=dat,
                                       contrasts=list("cut(x1, c(0, 3, 6))"=diag(2)))
    expect_equal(simple5@OriginalVariables, c(1,2,2,3))
    expect_equal(simple5@TermLabels, c("x2", "cut(x1, c(0, 3, 6))", "fac"))
    expect_equal(simple5@NM.Covariates, 1+c(0,1,1,2))
    expect_equal(simple5@NM.terms, 1+c(0,1,2))
 
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
    simple2 <- RItools:::design_matrix(z ~ x1 + I(x1^2) + fac, data=dat)
    expect_equal(simple2@OriginalVariables, 1:3)
    expect_equal(simple2@TermLabels, c("x1", "I(x1^2)", "fac"))
    expect_equal(simple2@NM.Covariates, 1+c(1,1,2))
    expect_equal(simple2@NM.terms, 1+c(1,1,2))
    ## If exactly two terms have missing data but in the same pattern, then
    ## NotMissing is a matrix w/ n rows and 1 col.
    simple3 <- RItools:::design_matrix(z ~ x1 + I(x1^2), data=dat)
    expect_equal(simple3@OriginalVariables, 1:2)
    expect_equal(simple3@TermLabels, c("x1", "I(x1^2)"))
    expect_equal(ncol(simple3@NotMissing), 2)
    expect_equal(simple3@NM.Covariates, 1 + c(1,1))
    expect_equal(simple3@NM.terms, 1 + c(1,1))

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

  # dropping the overall comparison
  expect_equal(dim(RItools:::makeDesigns(z.good ~ x + cluster(cluster) + strata(strata.good) - 1, data = d)@StrataFrame)[2], 1)
   
  ## More tests to write:
  # - All NA strata variables
  # - Missing z or cluster
  # - strata with extra levels but no observations (which can be safely dropped)
  #   (NB: extra levels tested upstream, in xBalance, as of commit 34861515; 
  #   see ./test.clusters.R ) 
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
  expect_equal(mean(d$x[d$z == 1]), descriptives["x", "Treatment", "Unstrat"])
  expect_equal(mean(d$x[d$z == 0]), descriptives["x", "Control", "Unstrat"])

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
  expect_equal(descriptives.all["x", "Treatment", "Unstrat"], mean(d$x[d$z == 1]))
  expect_equal(descriptives.all["x", "Treatment", "s"], mean(d$x[d$z == 1 & !is.na(d$s)]))

  simple.missing <- RItools:::makeDesigns(z ~ x + f + strata(s) + cluster(c), data = d.missing)
  
  descriptives.missing <- RItools:::designToDescriptives(simple.missing)

  with(d.missing,
       expect_equal(descriptives.missing["x", "Treatment", "Unstrat"],
                    mean(x[z == 1], na.rm = TRUE)))

  with(d.missing,
       expect_equal(descriptives.missing["x", "Treatment", "s"],
                    mean(x[z == 1 & !is.na(s)], na.rm = TRUE)))

  expect_false(identical(descriptives.all, descriptives.missing))

  # ETT weighting
  design.paired   <- RItools:::makeDesigns(z ~ x + f + strata(paired) + strata(s), data = d) 
  descriptives.paired <- RItools:::designToDescriptives(design.paired)

  with(d, expect_equal(descriptives.paired["x", "Control", "paired"], mean(x[z == 0])))
  with(d, expect_equal(descriptives.paired["x", "Control", "Unstrat"], mean(x[z == 0])))
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
    expect_equal(aggDesign@Covariates[,'x'], c(0, 0, 1, 1) )
    nm.column.for.x <- aggDesign@NM.Covariates[match( 'x', colnames(aggDesign@Covariates))]
  expect_equal(aggDesign@NotMissing[,nm.column.for.x], c(0,0,1,1) )
})


test_that("Aggregation of element weights to cluster level",{
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

  expect_equal(2*aggDesign.tall@ElementWeights,aggDesign.d2@ElementWeights)
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
     expect_equivalent(model.matrix(fff, m), as.matrix(design_matrix(fff, trees2)))
     trees2[1, "Height"] <- NA # RHS variable, but still shouldn't cause rows to be dropped
     expect_equal(dim(model.matrix(fff, m)), dim(as.matrix(design_matrix(fff, trees2))))


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

### Tests to write...
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
    dsimple2 <- RItools:::designToDescriptives(simple2)
    asimple2 <- RItools:::alignDesignsByStrata(simple2)
    expect_equivalent(dimnames(dsimple2)[[1]], colnames(asimple2[[1]]@Covariates))

})

##test_that("alignDesigns properly tracks ElementWeights vs NotMissing",{})
##test_that("",{})
