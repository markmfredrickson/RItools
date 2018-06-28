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
      expect_setequal(dimnames(xb[['results']])[['stat']],
                      c("Control", "Treatment","std.diff", "adj.diff", "pooled.sd", "z", "p")
                      ) #if any of these fail just have to adjust `tidy.xbal()`
      expect_equal(intersect(1:5, 4:2), 2:4) # order from 1st arg not 2nd
      expect_equal(dimnames(xb[['results']])[['strata']][1],'s') # "s" comes before, 
      expect_equal(dimnames(xb[['overall']])[[1]][1],'s') # not after, "--"
      expect_true(all( attr(xb[['results']], "NMpatterns") %in%
                       c("", dimnames(xb[['results']])[['vars']])
                      ) #usually all of the NM patterns should be names of vars,
                  )     #w/ exception of "", which indicates nothing is missing...
      ##...that is unless there are rows with all Xs being NA.
      xb2 <- balanceTest(z~x2,data=dat)
      expect_equal(attr(xb2[['results']], "NMpatterns"), c("(_non-null record_)", ""))
})

test_that("Basic function of xbal tidy() and glance() methods", {
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )
      dat[2,'x2'] <- NA

      xb1 <- balanceTest(z~x1, data=dat)
      expect_s3_class(tidy.xbal(xb1), 'data.frame')
      expect_equal(tidy.xbal(xb1)$`NA.info`, "")
      expect_setequal(colnames(tidy.xbal(xb1)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","NA.info", "statistic", "p.value"))
            expect_setequal(colnames(tidy.xbal(xb1, varnames_crosswalk=NULL)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","NA.info", "z", "p"))
                  expect_setequal(colnames(tidy.xbal(xb1, varnames_crosswalk=c('z'='Z', 'p'='P'))),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","NA.info", "Z", "P"))

      xb2 <- balanceTest(z~x1+x2+strata(s), data=dat)
      expect_equal(tidy.xbal(xb2)$`NA.info`, c("","(x2)", ""))
      
      expect_s3_class(glance.xbal(xb1), 'data.frame')
      expect_equal(rownames(glance.xbal(xb1, strata="--")), "--")
      expect_equal(rownames(glance.xbal(xb2, strata="--")), "--")
      expect_equal(rownames(glance.xbal(xb2, strata="s")), "s")

})
