################################################################################
# balanceTest() against a supplied covariance matrix (balanceSpec)
################################################################################
library("testthat")
context("balanceSpec: balance tests against a supplied covariance")

make_data <- function(n = 240, seed = 20261029) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- 0.6 * x1 + rnorm(n)
  x3 <- rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  s <- factor(rep(1:4, length.out = n))
  z <- rbinom(n, 1, 0.5)
  data.frame(z = z, x1 = x1, x2 = x2, x3 = x3, g = g, s = s,
             stringsAsFactors = FALSE)
}

dat <- make_data()
Sigma <- cov(expandCovariates(~ x1 + x2 + g, dat))          # singular: ga + gb + gc = 1
SigmaN <- cov(expandCovariates(~ x1 + x2 + x3, dat))        # full rank

test_that("expandCovariates gives the columns balanceTest tests", {
  x <- expandCovariates(~ x1 + x2 + g, dat)
  expect_equal(colnames(x), c("x1", "x2", "ga", "gb", "gc"))
  expect_equal(unname(x[, "x1"]), dat$x1)
  expect_equal(unname(x[, "ga"]), as.numeric(dat$g == "a"))
  bt <- balanceTest(z ~ x1 + x2 + g, data = dat)
  vars <- dimnames(bt$results)[["vars"]]
  expect_equal(colnames(x), vars[!startsWith(vars, "(")])
  expect_equal(attr(x, "originals"), c(1L, 2L, 3L, 3L, 3L))

  ## response, strata() and cluster() terms are ignored
  expect_equal(expandCovariates(z ~ x1 + g + strata(s), dat),
               expandCovariates(~ x1 + g, dat))
  expect_equal(expandCovariates(~ x1 + cluster(s) + g, dat),
               expandCovariates(~ x1 + g, dat))

  ## character columns become factors; NAs are kept
  d2 <- dat
  d2$g <- as.character(d2$g)
  d2$x1[3] <- NA
  x2 <- expandCovariates(~ x1 + g, d2)
  expect_equal(colnames(x2), c("x1", "ga", "gb", "gc"))
  expect_true(is.na(x2[3, "x1"]))
})

test_that("balanceTest is a generic and the formula interface is unchanged", {
  expect_true(is.function(balanceTest.formula))
  bt <- balanceTest(z ~ x1 + x2 + g + strata(s), data = dat)
  expect_s3_class(bt, "balancetest")
  expect_equal(bt, balanceTest(z ~ x1 + x2 + g + strata(s), dat))
  expect_equal(bt, balanceTest(data = dat, z ~ x1 + x2 + g + strata(s)))
  expect_equal(bt, balanceTest(fmla = z ~ x1 + x2 + g + strata(s), data = dat))
  ## an argument neither method knows is an error, not silently dropped
  expect_error(balanceTest(z ~ x1, data = dat, p.adjust.methods = "none"), "p.adjust.methods")
  expect_error(balanceTest(balanceSpec("z", Sigma), data = dat, subset = dat$x1 > 0), "subset")
  expect_error(balanceTest(balanceSpec("z", Sigma), data = dat, unit.weights = 1), "unit.weights")
  expect_equal(dimnames(bt$results)[["strata"]], c("s", "--"))
  ## named arguments still reach the formula method
  d <- dat
  d$w <- 1
  bt2 <- balanceTest(z ~ x1 + x2, data = d, p.adjust.method = "none", unit.weights = w)
  expect_s3_class(bt2, "balancetest")
  expect_equal(unlist(bt2$overall["--", ]),
               unlist(balanceTest(z ~ x1 + x2, data = d, p.adjust.method = "none")$overall["--", ]))
})

test_that("balanceSpec validates its covariance", {
  expect_error(balanceSpec("z", covariance = Sigma[, -1]), "square")
  m <- Sigma; dimnames(m) <- NULL
  expect_error(balanceSpec("z", covariance = m), "names")
  m <- Sigma; m[1, 2] <- m[1, 2] + 1
  expect_error(balanceSpec("z", covariance = m), "symmetric")
  m <- Sigma; m[1, 1] <- NA
  expect_error(balanceSpec("z", covariance = m), "missing")
  m <- Sigma; m["x1", "x2"] <- m["x2", "x1"] <- 5 * sqrt(m["x1", "x1"] * m["x2", "x2"])
  expect_error(balanceSpec("z", covariance = m), "positive semi-definite")
  expect_error(balanceSpec("z", covariance = Sigma, inverse_sqrt = Sigma), "'r' attribute")
  inv <- XtX_pseudoinv_sqrt(Sigma, mat.is.XtX = TRUE)
  wrong <- inv[rev(seq_len(nrow(inv))), , drop = FALSE]
  rownames(wrong) <- rev(colnames(Sigma))
  attr(wrong, "r") <- attr(inv, "r")
  expect_error(balanceSpec("z", covariance = Sigma, inverse_sqrt = wrong), "order")
  expect_error(balanceTest(balanceSpec("z", Sigma), data = dat, unstratified = FALSE),
               "needs strata")
  spec <- balanceSpec("z", covariance = Sigma)
  expect_s3_class(spec, "balanceSpec")
  expect_equal(attr(spec$inverse_sqrt, "r"), 4L) # three indicator columns sum to one
  expect_output(print(spec), "rank 4")
})

test_that("with the sample covariance, the unstratified test reproduces HB08 exactly", {
  ## full-rank covariates: the omnibus statistic and its df must agree exactly
  spec <- balanceSpec(treatment = "z", covariance = SigmaN,
                      variables = c("x1", "x2", "x3"))
  bs <- balanceTest(spec, data = dat)
  bt <- balanceTest(z ~ x1 + x2 + x3, data = dat)
  expect_s3_class(bs, "balancetest")
  expect_equal(attr(bs, "covariance"), "supplied")
  expect_equal(dimnames(bs$results)[["strata"]], "--")
  expect_equal(bs$overall["--", ], bt$overall["--", ])
  expect_equal(bs$overall["--", "df"], 3)
  expect_equal(attr(bs$overall, "tcov")[["--"]], attr(bt$overall, "tcov")[["--"]])
  for (stat in c("adj.diff", "z", "p", "Treatment", "Control")) {
    expect_equal(bs$results[, stat, "--"], bt$results[, stat, "--"], info = stat)
  }
  ## standard differences are scaled by the supplied covariance's diagonal
  expect_equal(bs$results[c("x1", "x2"), "pooled.sd", "--"], sqrt(diag(SigmaN))[c("x1", "x2")])
  expect_equal(bs$results["x1", "std.diff", "--"],
               bs$results["x1", "adj.diff", "--"] / sqrt(SigmaN["x1", "x1"]))

  ## with a factor the covariance is singular (its indicators sum to one). The
  ## univariate statistics agree exactly, and the omnibus df is the rank, 4,
  ## on both paths: each gives XtX_pseudoinv_sqrt() a tolerance that cuts on
  ## the matrix's own scale, not at machine precision (#146)
  spec <- balanceSpec(treatment = "z", covariance = Sigma, variables = c("x1", "x2", "g"))
  bs <- balanceTest(spec, data = dat)
  bt <- balanceTest(z ~ x1 + x2 + g, data = dat)
  expect_equal(attr(bs$overall, "tcov")[["--"]], attr(bt$overall, "tcov")[["--"]])
  for (stat in c("adj.diff", "z", "p")) {
    expect_equal(bs$results[, stat, "--"], bt$results[, stat, "--"], info = stat)
  }
  expect_equal(bs$overall["--", "df"], 4)
  expect_equal(bt$overall["--", "df"], 4)
  expect_equal(bs$overall["--", "chisquare"], bt$overall["--", "chisquare"])
  expect_equal(attr(spec$inverse_sqrt, "r"), 4L)
})

test_that("the rank cut depends on neither rounding nor the variables' units (#146)", {
  ## rounding-level noise on a singular covariance does not raise its rank
  noisy <- Sigma + 1e-13 * tcrossprod(rep(1, nrow(Sigma)))
  expect_equal(attr(RItools:::.pseudoinvSqrtOfCovariance(noisy), "r"), 4)
  ## the formula method's omnibus df with a factor is the rank on both rows
  bt <- balanceTest(z ~ x1 + x2 + g + strata(s), data = dat)
  expect_equal(bt$overall[, "df"], c(4, 4))

  ## a count in the thousands beside a proportion: both directions are kept
  ## and the statistic does not depend on the units
  set.seed(146)
  d <- dat
  d$enroll <- round(1000 + 800 * d$x1 + 300 * rnorm(nrow(d)))
  d$prop <- plogis(d$x2) / 20
  bt <- balanceTest(z ~ enroll + prop, data = d)
  expect_equal(bt$overall["--", "df"], 2)
  d$enroll_k <- d$enroll / 1e6
  bt_k <- balanceTest(z ~ enroll_k + prop, data = d)
  expect_equal(bt_k$overall["--", "chisquare"], bt$overall["--", "chisquare"])
  S <- cov(expandCovariates(~ enroll + prop, d))
  expect_gt(S["enroll", "enroll"] / S["prop", "prop"], 1e8)
  spec <- balanceSpec("z", S, variables = c("enroll", "prop"))
  expect_equal(attr(spec$inverse_sqrt, "r"), 2)
  bs <- balanceTest(spec, data = d)
  expect_equal(bs$overall["--", "df"], 2)
  expect_equal(bs$overall["--", "chisquare"], bt$overall["--", "chisquare"])
  ## the wrapper returns a generalized inverse: V V' S V V' = V V'
  V <- RItools:::.pseudoinvSqrtOfCovariance(Sigma)
  G <- tcrossprod(V)
  expect_equal(Sigma %*% G %*% Sigma, Sigma, check.attributes = FALSE)
  ## a slightly negative eigenvalue from rounding is not counted as a direction
  e <- eigen(cov2cor(SigmaN)); lam <- e$values; lam[3] <- -1e-10
  nearly <- e$vectors %*% (lam * t(e$vectors))
  dimnames(nearly) <- dimnames(SigmaN)
  expect_equal(attr(balanceSpec("z", nearly)$inverse_sqrt, "r"), 2)
  bad <- SigmaN; bad[1, 1] <- -1
  expect_error(balanceSpec("z", bad), "negative variance")
})

test_that("the formula method's model frame holds the formula's columns (#146)", {
  ## unit.weights and subset expressions built from the caller's variables
  ## work inside a function, where they used to be looked up in the global
  ## environment only
  f <- function(d) {
    w <- rep(2, nrow(d))
    keep <- d$x1 > -10
    list(balanceTest(z ~ x1 + x2 + strata(s), data = d, unit.weights = w),
         balanceTest(z ~ x1 + x2 + strata(s), data = d, subset = keep),
         balanceTest(z ~ x1 + x2 + strata(s), data = d, unit.weights = rep(1, nrow(d))))
  }
  res <- f(dat)
  ref <- balanceTest(z ~ x1 + x2 + strata(s), data = dat)
  for (r in res) {
    expect_equal(r$overall[, "chisquare"], ref$overall[, "chisquare"])
    expect_equal(r$results[, "adj.diff", ], ref$results[, "adj.diff", ])
  }
  ## a subset given by a local vector takes effect
  g <- function(d) { keep <- d$g != "c"; balanceTest(z ~ x1 + x2, data = d, subset = keep) }
  expect_equal(g(dat)$results[, "adj.diff", "--"],
               balanceTest(z ~ x1 + x2, data = dat[dat$g != "c", ])$results[, "adj.diff", "--"])
  ## a term that uses an object that is not a column, of another length
  brks <- c(-Inf, 0, Inf)
  bt <- balanceTest(z ~ cut(x1, brks) + x2, data = dat)
  expect_equal(sum(grepl("^cut", dimnames(bt$results)[["vars"]])), 2)
  ## `.` still means every column, and a column the formula does not name
  ## does not matter
  d2 <- dat[, c("z", "x1", "x2")]
  expect_equal(balanceTest(z ~ ., data = d2)$overall, balanceTest(z ~ x1 + x2, data = d2)$overall)
  d3 <- dat; d3$junk <- list(1)[rep(1, nrow(d3))]
  expect_equal(balanceTest(z ~ x1 + x2, data = d3)$overall, balanceTest(z ~ x1 + x2, data = dat)$overall)
  ## a column whose name is not a syntactic name
  d4 <- dat; d4$`x 3` <- d4$x3; d4$"a`b" <- d4$x2
  expect_equal(unname(balanceTest(z ~ `x 3`, data = d4)$results[1, "adj.diff", "--"]),
               unname(balanceTest(z ~ x3, data = dat)$results[1, "adj.diff", "--"]))
  expect_s3_class(balanceTest(z ~ ., data = d4[, c("z", "x1", "a`b")]), "balancetest")
  ## a formula with no environment still works
  fm <- z ~ x1 + x2; environment(fm) <- NULL
  expect_equal(balanceTest(fm, data = dat)$overall, balanceTest(z ~ x1 + x2, data = dat)$overall)
})

test_that("stratified: exact agreement with HB08 when every stratum has the same covariance", {
  ## the same covariate block in each of four strata, shifted, so the
  ## within-stratum sample covariances coincide
  set.seed(1)
  block <- cbind(x1 = rnorm(60), x2 = rnorm(60))
  block[, 2] <- 0.5 * block[, 1] + block[, 2]
  d <- do.call(rbind, lapply(1:4, function(b) {
    data.frame(x1 = block[, 1] + b, x2 = block[, 2] - b, s = b)
  }))
  d$s <- factor(d$s)
  repeat {
    d$z <- rbinom(nrow(d), 1, 0.4)
    if (all(tapply(d$z, d$s, function(v) length(unique(v)) == 2))) break
  }
  S <- cov(block)
  spec <- balanceSpec("z", covariance = S, strata = "s", variables = c("x1", "x2"))
  bs <- balanceTest(spec, data = d)
  bt <- balanceTest(z ~ x1 + x2 + strata(s), data = d)
  expect_equal(dimnames(bs$results)[["strata"]], c("s", "--"))
  expect_equal(unlist(bs$overall["s", ]), unlist(bt$overall["s", ]))
  expect_equal(attr(bs$overall, "tcov")[["s"]], attr(bt$overall, "tcov")[["s"]])
  expect_equal(bs$results[, "z", "s"], bt$results[, "z", "s"])
  expect_equal(bs$results[, "adj.diff", "s"], bt$results[, "adj.diff", "s"])
  ## the unstratified row is not exact here (stratum means differ) but is present
  expect_true(is.finite(bs$overall["--", "chisquare"]))

  ## unstratified = FALSE drops it
  bs1 <- balanceTest(spec, data = d, unstratified = FALSE)
  expect_equal(dimnames(bs1$results)[["strata"]], "s")
  expect_equal(unlist(bs1$overall["s", ]), unlist(bs$overall["s", ]))

  ## vectors in place of column names give the same answer
  spec_v <- balanceSpec(d$z, covariance = S, strata = d$s, variables = c("x1", "x2"))
  bs_v <- balanceTest(spec_v, data = d[, c("x1", "x2")])
  expect_equal(unname(unlist(bs_v$overall[1, ])), unname(unlist(bs$overall["s", ])))
  expect_equal(dimnames(bs_v$results)[["strata"]], c(".strata", "--"))
})

test_that("with the true covariance, the stratified test is close to HB08 on a large sample", {
  set.seed(2)
  n <- 4000
  s <- factor(rep(1:8, length.out = n))
  x1 <- rnorm(n) + as.numeric(s)
  x2 <- 0.7 * x1 + rnorm(n)
  z <- rbinom(n, 1, 0.5)
  d <- data.frame(z, x1, x2, s)
  Strue <- matrix(c(1, 0.7, 0.7, 0.49 + 1), 2, 2, dimnames = list(c("x1", "x2"), c("x1", "x2")))
  bs <- balanceTest(balanceSpec("z", Strue, strata = "s"), data = d)
  bt <- balanceTest(z ~ x1 + x2 + strata(s), data = d)
  ratio <- attr(bs$overall, "tcov")[["s"]] / attr(bt$overall, "tcov")[["s"]]
  expect_true(all(abs(ratio - 1) < 0.1))
  expect_equal(bs$overall["s", "df"], bt$overall["s", "df"])
})

test_that("variables default to every column not otherwise used", {
  d <- dat[, setdiff(names(dat), "x3")]
  spec <- balanceSpec("z", covariance = Sigma, strata = "s")
  bs <- balanceTest(spec, data = d)
  vars <- dimnames(bs$results)[["vars"]]
  expect_equal(vars[!startsWith(vars, "(")], colnames(Sigma))
  ## and a one-sided formula works
  spec_f <- balanceSpec("z", covariance = Sigma, strata = "s", variables = ~ x1 + x2 + g)
  expect_equal(balanceTest(spec_f, data = d)$overall, bs$overall)
  ## a column the covariance lacks is an error, not silently untested
  expect_error(balanceTest(spec, data = dat), "x3")
})

test_that("an expanded covariate column missing from the covariance is an error", {
  spec <- balanceSpec("z", covariance = Sigma[c("x1", "x2"), c("x1", "x2")],
                      variables = c("x1", "x2", "g"))
  expect_error(balanceTest(spec, data = dat), "ga, gb, gc")
  spec2 <- balanceSpec("z", covariance = Sigma, variables = c("x1", "nope"))
  expect_error(balanceTest(spec2, data = dat), "nope")
  expect_error(balanceTest(balanceSpec("zz", covariance = Sigma), data = dat), "treatment")
})

test_that("covariance columns absent from the data are set aside with a message", {
  d <- dat[dat$g != "c", ]
  spec <- balanceSpec("z", covariance = Sigma, variables = c("x1", "x2", "g"), strata = "s")
  msgs <- capture_messages(bs <- balanceTest(spec, data = d))
  expect_length(msgs, 1L)   # once, not once per stratification
  expect_match(msgs, "gc")
  vars <- dimnames(bs$results)[["vars"]]
  expect_false("gc" %in% vars)
  S2 <- Sigma[c("x1", "x2", "ga", "gb"), c("x1", "x2", "ga", "gb")]
  bs2 <- balanceTest(balanceSpec("z", covariance = S2, variables = c("x1", "x2", "g"), strata = "s"),
                     data = d)
  expect_equal(bs$overall, bs2$overall)
})

test_that("a zero entry on the covariance diagonal leaves that scale to the sample", {
  m <- Sigma
  m["x2", ] <- m[, "x2"] <- 0
  bs <- balanceTest(balanceSpec("z", covariance = m, variables = c("x1", "x2", "g")), data = dat)
  bt <- balanceTest(z ~ x1 + x2 + g, data = dat)
  expect_equal(bs$results["x2", "pooled.sd", "--"], bt$results["x2", "pooled.sd", "--"])
  expect_true(is.finite(bs$results["x2", "std.diff", "--"]))
  expect_true(is.na(bs$results["x2", "z", "--"]))        # no null variance supplied for it
  expect_equal(bs$overall["--", "df"], 3)                  # x1 plus two for the factor
})

test_that("a precomputed inverse square root is used as given", {
  ## (full-rank matrix: XtX_pseudoinv_sqrt()'s own rank cut sits at machine
  ## precision and need not agree with balanceSpec()'s on a singular one)
  inv <- XtX_pseudoinv_sqrt(SigmaN, mat.is.XtX = TRUE)
  vars <- c("x1", "x2", "x3")
  spec <- balanceSpec("z", covariance = SigmaN, variables = vars, inverse_sqrt = inv)
  spec0 <- balanceSpec("z", covariance = SigmaN, variables = vars)
  expect_equal(balanceTest(spec, data = dat), balanceTest(spec0, data = dat))
  ## a wrong inverse changes the omnibus statistic but nothing else
  bad <- inv * 2
  attr(bad, "r") <- attr(inv, "r")
  bs_bad <- balanceTest(balanceSpec("z", covariance = SigmaN, variables = vars,
                                    inverse_sqrt = bad), data = dat)
  bs <- balanceTest(spec0, data = dat)
  expect_equal(bs_bad$overall["--", "chisquare"], 4 * bs$overall["--", "chisquare"])
  expect_equal(bs_bad$results[, "z", "--"], bs$results[, "z", "--"])
})

test_that("clusters and unit weights pass through to the design", {
  d <- dat
  d$s <- factor(rep(1:4, each = nrow(d) / 4))          # strata of 60 rows
  d$cl <- factor(rep(seq_len(nrow(d) / 4), each = 4))  # clusters of 4 rows, nested in s
  d$z <- as.numeric(as.integer(d$cl) %% 2 == 0)
  d$w <- runif(nrow(d), 0.5, 2)
  spec <- balanceSpec("z", covariance = Sigma[c("x1", "x2"), c("x1", "x2")],
                      strata = "s", cluster = "cl", unit.weights = "w",
                      variables = c("x1", "x2"))
  bs <- balanceTest(spec, data = d)
  expect_equal(dimnames(bs$results)[["strata"]], c("s", "--"))
  expect_true(all(is.finite(bs$overall[, "chisquare"])))
  expect_warning(balanceTest(balanceSpec("z", covariance = Sigma[1:2, 1:2], strata = "s",
                                         unit.weights = ifelse(seq_len(nrow(d)) == 1, NA, 1),
                                         variables = c("x1", "x2")), data = d),
                 "NA unit.weights")
})

test_that("missing covariate values reduce the null variance as in the formula method", {
  d <- dat
  d$x1[1:10] <- NA
  spec <- balanceSpec("z", covariance = Sigma, variables = c("x1", "x2", "g"))
  bs <- balanceTest(spec, data = d)
  vars <- dimnames(bs$results)[["vars"]]
  expect_true(any(startsWith(vars, "(")))
  expect_true(all(is.na(bs$results[startsWith(vars, "("), "z", "--"])))
  expect_true(is.finite(bs$overall["--", "chisquare"]))

  ## a small sample with as many tested columns as units warns, as HB08 does
  tiny <- dat[1:5, ]          # rank 4 against 5 units less one stratum
  tiny$z <- c(1, 0, 1, 0, 1)
  expect_warning(balanceTest(balanceSpec("z", Sigma, variables = c("x1", "x2", "g")), data = tiny),
                 "degenerate")

  ## with 30% of x1 missing at random on a large sample, the supplied-covariance
  ## null variance of x1 tracks HB08's, and is about 0.7 of the complete-data one
  set.seed(3)
  n <- 6000
  big <- data.frame(z = rbinom(n, 1, 0.5), x1 = rnorm(n), s = factor(rep(1:6, length.out = n)))
  big$x2 <- 0.6 * big$x1 + rnorm(n)
  Strue <- matrix(c(1, 0.6, 0.6, 1.36), 2, 2, dimnames = list(c("x1", "x2"), c("x1", "x2")))
  full <- balanceTest(balanceSpec("z", Strue, strata = "s", variables = c("x1", "x2")), data = big)
  big$x1[sample(n, 0.3 * n)] <- NA
  miss <- balanceTest(balanceSpec("z", Strue, strata = "s", variables = c("x1", "x2")), data = big)
  ref <- balanceTest(z ~ x1 + x2 + strata(s), data = big)
  tm <- attr(miss$overall, "tcov")[["s"]][c("x1", "x2"), c("x1", "x2")]
  tr <- attr(ref$overall, "tcov")[["s"]][c("x1", "x2"), c("x1", "x2")]
  tf <- attr(full$overall, "tcov")[["s"]][c("x1", "x2"), c("x1", "x2")]
  expect_true(all(abs(tm / tr - 1) < 0.1))
  expect_equal(tm["x1", "x1"] / tf["x1", "x1"], 0.7, tolerance = 0.03)
  expect_equal(tm["x2", "x2"], tf["x2", "x2"])
  expect_equal(miss$overall["s", "df"], 2)
})
