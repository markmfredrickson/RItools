library("testthat")

context("balanceTest plots")

test_that("Uses ggplot", {

  set.seed(20121119)
  Z <- rep(c(0,1), 10)
  df <- data.frame(Z = Z,
                   X1 = rnorm(20, mean = Z*2),
                   X2 = rnorm(20, mean = Z*3),
                   X3 = rnorm(20, mean = Z * -1),
                   X4 = rnorm(20, mean = Z * -0.5),
                   s = rep(1:4, each = 5))

  bt <- balanceTest(Z ~ X1 + X2 + X3 + X4 + strata(s), data = df)
  btp <- plot(bt)
  expect_is(btp, "gg")
  btp

})


if (requireNamespace("optmatch", quietly = TRUE)) {
  test_that("plotting with nuclearplants", {
    data(nuclearplants)
    psm <- glm(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + pt,
           family = binomial, data = nuclearplants)
    psm.dist <- optmatch::match_on(psm, data=nuclearplants)
    ps.pm2 <- optmatch::pairmatch(psm.dist, data = nuclearplants)
    myb <- balanceTest(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n +
                         strata(ps.pm2), data = nuclearplants)
    p <- plot(myb)
    expect_is(p, "gg")
    p

  })
}
