context("broom package-style model tidiers")

test_that("tidy.xbal assumptions haven't changed",{
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )

      xb <- balanceTest(z~x1+strata(s), data=dat)
      expect_equal(length(dim(xb[['results']])), 3)
      expect_equal(names(dimnames(xb[['results']]))[1:2],
                   c("vars", "stat") )
      expect_setequal(dimnames(xb[['results']])[['stat']],
                      c("Control", "Treatment","std.diff", "adj.diff", "pooled.sd", "z", "p")
                      ) #if any of these fail just have to adjust `tidy.xbal()`
      expect_equal(intersect(1:5, 4:2), 2:4) # order from 1st arg not 2nd
})

test_that("Basic function of xbal tidy() and glance() methods", {
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )

      xb <- balanceTest(z~x1+strata(s), data=dat)
      expect_s3_class(tidy.xbal(xb), 'data.frame')
      expect_setequal(colnames(tidy.xbal(xb)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","statistic", "p.value"))
            expect_setequal(colnames(tidy.xbal(xb, varnames_crosswalk=NULL)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","z", "p"))
                  expect_setequal(colnames(tidy.xbal(xb, varnames_crosswalk=c('z'='Z', 'p'='P'))),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","Z", "P"))
            expect_setequal(colnames(tidy.xbal(xb)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","statistic", "p.value"))

      expect_s3_class(glance.xbal(xb), 'data.frame')
      expect_equal(rownames(glance.xbal(xb, strata="s")), "s")
            expect_equal(rownames(glance.xbal(xb, strata="Unstrat")), "Unstrat")
})
