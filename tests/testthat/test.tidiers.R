context("broom package-style model tidiers")

test_that("tidy.xbal assumptions haven't changed",{
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )
      dat[2,'x2'] <- NA

      xb <- balanceTest(z~x1+x2+strata(s), data=dat)
      expect_equal(length(dim(xb[['results']])), 3)
      expect_equal(names(dimnames(xb[['results']]))[1:2],
                   c("vars", "stat") )
      expect_true(setequal(dimnames(xb[['results']])[['stat']],
                      c("Control", "Treatment","std.diff", "adj.diff", "pooled.sd", "z", "p")
                      )) #if any of these fail just have to adjust `tidy.xbal()`
      expect_equal(intersect(1:5, 4:2), 2:4) # order from 1st arg not 2nd
      expect_equal(dimnames(xb[['results']])[['strata']][1],'s') # "s" comes before, 
      expect_equal(dimnames(xb[['overall']])[[1]][1],'s') # not after, "--"
      expect_true(all( attr(xb[['results']], "NMpatterns") %in%
                       c("", dimnames(xb[['results']])[['vars']])
                      ) #usually all of the NM patterns should be names of vars,
                  )     #w/ exception of "", which indicates nothing is missing...
      ##...that is unless there are rows with all Xs being NA.
      xb2 <- balanceTest(z~x2,data=dat)
      expect_equal(attr(xb2[['results']], "NMpatterns"), c("(_any Xs recorded_)", ""))
})

test_that("Basic function of xbal tidy() and glance() methods", {
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )
      dat[2,'x2'] <- NA

      xb <- balanceTest(z~x1+strata(s), data=dat)
      expect_s3_class(tidy.xbal(xb), 'data.frame')
      expect_true(setequal(colnames(tidy.xbal(xb)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd", "NA.info", "statistic", "p.value")))
            expect_true(setequal(colnames(tidy.xbal(xb, varnames_crosswalk=NULL)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","z", "p", "NA.info")))
                  expect_true(setequal(colnames(tidy.xbal(xb, varnames_crosswalk=c('z'='Z', 'p'='P'))),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","Z", "P", "NA.info")))

})

test_that("tidy.xbal w/ special formatting for original units vars",{
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )
      dat[2,'x2'] <- NA

      xb1 <- balanceTest(z~x1, data=dat)
      t1 <- tidy.xbal(xb1, format=TRUE, digits=2)
      expect_s3_class(t1, 'data.frame')
      expect_is(t1$"adj.diff", "character")
    })
          
