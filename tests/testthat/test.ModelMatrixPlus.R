context("ModelMatrixPlus")



test_that("Model matrix material is properly formed", {

     ff <- log(Volume) ~ log(Height) + log(Girth)
     trees1 <- trees
     trees1$'(weights)' <- 1
     DM0 <- model_matrix(ff, trees1, remove.intercept=FALSE)
     expect_is(DM0, "ModelMatrixPlus")
     m <- model.frame(ff, trees)
     expect_equivalent(model.matrix(ff, m), as.matrix(DM0))
     fff <- update(ff, .~.-1)
     trees2 <- trees
     trees2[1, "Volume"] <- NA # LHS variable, shouldn't affect anything
     m2a <- model.frame(fff, trees2, na.action = na.pass)
     m2a$'(weights)' <- 1
     expect_equivalent(model.matrix(fff, m), as.matrix(model_matrix(fff, m2a)))
     trees2[1, "Height"] <- NA # RHS variable, but still shouldn't cause rows to be dropped
     m2b <- model.frame(fff, trees2, na.action = na.pass)
     m2b$'(weights)' <- 1
     expect_equal(dim(model.matrix(fff, m)), 
                  dim(as.matrix(model_matrix(fff,m2b) ) )
                  )

     ## specified contrasts
     dd <- data.frame(a = gl(3,4), b = gl(4,1,12)) # balanced 2-way
     dd$'(weights)' <- 1
     expect_equal(model.matrix(~ a + b, dd),
                  as.matrix(model_matrix(~ a + b, dd, remove.intercept=FALSE)))
     expect_equal(model.matrix(~ a + b-1, dd),
                  as.matrix(model_matrix(~ a + b-1, dd)))
     expect_equal(model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum")),
                  as.matrix(model_matrix(~ a + b, dd,
                                          contrasts = list(a = "contr.sum"), remove.intercept=FALSE))
                  )
     expect_equal(model.matrix(~ a + b, dd,
                               contrasts = list(a = "contr.sum", b = "contr.poly")),
                  as.matrix(model_matrix(~ a + b, dd,
                                          contrasts =
                                              list(a = "contr.sum", b = "contr.poly"),
                                          remove.intercept=FALSE))
                  )
          }
          )

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
  simple <- RItools:::model_matrix(z.good ~ x, data = d)
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
              datmf <- model.frame(z ~ x1 + x2 + fac, dat, na.action = na.pass) 
              datmf$'(weights)' <- 1
              simple2 <- RItools:::model_matrix(z ~ x1 + x2 + fac, data = datmf)
              expect_equivalent(ncol(simple2@NotMissing), 3)
              expect_equivalent(colnames(simple2@NotMissing), c("_non-null record_", "x1", "fac"))
              
})
test_that("lookup tables OK, even w/ complex & multi-column terms",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=rep(c(NA, 1), c(3,2)),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
    datmf <- model.frame(z ~ x1 + x2 + fac, dat, na.action = na.pass)
    datmf$'(weights)' <- 1
    simple2 <- RItools:::model_matrix(z ~ x1 + x2 + fac, data=datmf)
    expect_equal(simple2@OriginalVariables, 1:3)
    expect_equal(simple2@TermLabels, c( "x1", "x2", "fac"))
    expect_equal(simple2@NM.Covariates, c(2,0,3))
    expect_equal(simple2@NM.terms, c(2,0,3))

    ## check that complex term don't spell trouble in themselves
    datmf <- model.frame(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = dat,
                         na.action = na.pass)
    datmf$'(weights)' <- 1    
    simple3 <- RItools:::model_matrix(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = datmf)
    expect_equal(simple3@OriginalVariables, 1:3)
    expect_equal(simple3@TermLabels, c("x1", "cut(x2, c(0, 3, 6))", "fac"))
    expect_equal(simple3@NM.Covariates, c(2,0,3))
    expect_equal(simple3@NM.terms, c(2,0,3))

    ## now try a complex term that actually expands to multiple columns
    datmf <- model.frame(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = dat,
                         na.action = na.pass)
    datmf$'(weights)' <- 1
    simple4 <- RItools:::model_matrix(z ~ x1 + cut(x2, c(0,3,6)) + fac, data = datmf,
                                       contrasts=list("cut(x2, c(0, 3, 6))"=diag(2)))
    expect_equal(simple4@OriginalVariables, c(1,2,2,3))
    expect_equal(simple4@TermLabels, c("x1", "cut(x2, c(0, 3, 6))", "fac"))
    expect_equal(simple4@NM.Covariates, c(2,0,0,3))
    expect_equal(simple4@NM.terms, c(2,0,3))

    ## now put NAs in the multi-column complex term
    datmf <- model.frame(z ~ x2 + cut(x1, c(0,3,6)) + fac, data = dat,
                         na.action = na.pass)
    datmf$'(weights)' <- 1
    simple5 <- RItools:::model_matrix(z ~ x2 + cut(x1, c(0,3,6)) + fac, data = datmf,
                                       contrasts=list("cut(x1, c(0, 3, 6))"=diag(2)))
    expect_equal(simple5@OriginalVariables, c(1,2,2,3))
    expect_equal(simple5@TermLabels, c("x2", "cut(x1, c(0, 3, 6))", "fac"))
    expect_equal(simple5@NM.Covariates, c(0,2,2,3))
    expect_equal(simple5@NM.terms, c(0,2,3))
 
})


test_that("Duplicated missingness patterns handled appropriately",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=c(rep(NA, 3), 1:2),
                    x2=c(1:5),
                    fac=factor(c(rep(1:2,2), NA))
                    )
    datmf <- model.frame(z ~ x1 + I(x1^2) + fac, data=dat, na.action = na.pass)
    datmf$'(weights)' <- 1
    simple2 <- RItools:::model_matrix(z ~ x1 + I(x1^2) + fac, data = datmf)
    expect_equal(simple2@OriginalVariables, 1:3)
    expect_equal(simple2@TermLabels, c("x1", "I(x1^2)", "fac"))
    expect_equal(simple2@NM.Covariates, 1+c(1,1,2))
    expect_equal(simple2@NM.terms, 1+c(1,1,2))

    ## If exactly two terms have missing data but in the same pattern, then
    ## NotMissing is a matrix w/ n rows and 1 col.
    datmf <- model.frame(z ~ x1 + I(x1^2), data=dat, na.action = na.pass)
    datmf$'(weights)' <- 1
    simple3 <- RItools:::model_matrix(z ~ x1 + I(x1^2), data = datmf)
    expect_equal(simple3@OriginalVariables, 1:2)
    expect_equal(simple3@TermLabels, c("x1", "I(x1^2)"))
    expect_equal(ncol(simple3@NotMissing), 1)
    expect_equal(simple3@NM.Covariates, c(1,1))
    expect_equal(simple3@NM.terms, c(1,1))

})

test_that("All-fields missingness |-> NotMissing col '_non-null record_'",{

    dat <- data.frame(strat=rep(letters[1:2], c(3,2)),
                    clus=factor(c(1,1,2:4)),
                    z=c(TRUE, rep(c(TRUE, FALSE), 2)),
                    x1=c(rep(NA, 4), 2),
                    x2=c(1:3, NA, 5),
                    fac=factor(c(NA, 1:2, NA, 2))
                    )
    dat$'(weights)' <- 1
    datmf <- model.frame(z ~ x1 + I(x1^2) + fac, data=dat, na.action = na.pass)
    datmf$'(weights)' <- 1
    simple6 <- RItools:::model_matrix(z ~ x1 + I(x1^2) + fac, data = datmf)
    expect_equal(simple6@OriginalVariables, 1:3)
    expect_equal(simple6@TermLabels, c("x1", "I(x1^2)", "fac"))
    expect_equal(colnames(simple6@NotMissing)[1], "_non-null record_")
    expect_equal(simple6@NM.Covariates, c(2,2,1))
    expect_equal(simple6@NM.terms, c(2,2,1))


} )
