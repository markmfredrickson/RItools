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
                       c("(_any Xs recorded_)", "", dimnames(xb[['results']])[['vars']])
                      ) #usually all of the NM patterns should be names of vars, w/ exceptions
                  )     #of "(_any Xs recorded_)" and  "", which indicates nothing is missing.
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

      bt <- balanceTest(z~x1 + strata(s), data=dat)
      expect_s3_class(tidy.xbal(bt), 'data.frame')
      expect_true(setequal(colnames(tidy.xbal(bt)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd", "NA.info", "statistic", "p.value")))
      expect_true(setequal(colnames(tidy.xbal(bt, varnames_crosswalk=NULL)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","z", "p", "NA.info")))
      expect_true(setequal(colnames(tidy.xbal(bt, varnames_crosswalk=c('z'='Z', 'p'='P'))),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","Z", "P", "NA.info")))
      expect_true(setequal(colnames(tidy.xbal(bt, varnames_crosswalk=c('Control'='notZ', 'p'='P'))),
                      c("vars","notZ", "Treatment","std.diff",
                        "adj.diff", "pooled.sd","z", "P", "NA.info"))
                  )
      expect_equal(colnames(tidy.xbal(bt, varnames_crosswalk=c('Control'='notZ', 'p'='P'))),
                   c("vars","notZ", "Treatment","std.diff",
                     "adj.diff", "pooled.sd","z", "P", "NA.info")
                   )
})

test_that("basic functionality of tidy.xbal as applied to value of a call to xBal()",{
      set.seed(20160406)
      n <- 7 
      dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
      dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )
      dat[2,'x2'] <- NA

      xb <- xBalance(z~x1, report="all", data=dat)
      expect_false(inherits(try(tidy.xbal(xb)), "try-error"))

      expect_s3_class(tidy.xbal(xb), 'data.frame')
      expect_true(setequal(colnames(tidy.xbal(xb)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff", "NA.info", "statistic", "p.value")))

      expect_true(setequal(colnames(tidy.xbal(xb, varnames_crosswalk=NULL)),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff",  "NA.info", "z", "p")))
                  expect_true(setequal(colnames(tidy.xbal(xb, varnames_crosswalk=c('z'='Z', 'p'='P'))),
                      c("vars","Control", "Treatment","std.diff",
                        "adj.diff",  "NA.info", "Z", "P")))

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
          

test_that("Date presentation", {
    set.seed(39483293)

    n <- 100
    z = rbinom(n, size = 1, p = 0.3)
    mydata <- data.frame(
        z = z,
        x1 = rnorm(n),
        x2 = as.Date(
            ifelse(z == 1, as.Date("2018-01-01"), as.Date("2017-01-01"))
            + rpois(n = n, lambda = 100), origin = "1970-01-01")
    )

    bt <- balanceTest(z ~ x1 + x2, data = mydata)

    library(broom)

    btt <- tidy(bt, format = TRUE,
                var_format = list("x2" = list(mean = function(x) { as.Date(x, origin = "1970-01-01") },
                                              diff = function(x) { as.difftime(x, units = "days")}))
                                       )

    expect_equal(as.character(btt[2, "Treatment"]), "2018-04-08")
    expect_equal(as.character(btt[2, "adj.diff"]), "366 days")

})
