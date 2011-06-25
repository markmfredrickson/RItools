################################################################################
# Fisher's Exact test, with classic examples
################################################################################

library(testthat)
context("Fisher exact test")

test_that("Lady Tasting Tea", {
  actual.cups <- c(1, 0, 0, 1, 1, 0, 1, 0)
  lady.guess <- c(1, 0, 0, 1, 1, 0, 0, 1)

  number.correct <- function(guesses, z, blocks) { sum(z == guesses) / 2 }
  expect_equal(number.correct(lady.guess, actual.cups, NULL), 3)
  
  lady.distribution <- parameterizedRandomizationDistribution(lady.guess, actual.cups,
    test.stat = number.correct)

  pv <- upper.p.value(lady.distribution[1,1],
    lady.distribution[1,-1])
  expect_equal(pv, 17/70)

  expect_equivalent(as.numeric(table(lady.distribution[1,-1])/70), 
    c(1/70, 16/70, 36/70, 16/70, 1/70))
  
})

test_that("2x2 table style", {
  set.seed(20110620)
  tau <- 10
  n <- 12
  Yc <- rbinom(p = 0.25, size = 1, n = n) 
  Yt <- Yc  + rbinom(p = 0.75, size = 1, n = n)
  Yt <- sapply(Yt, function(i) { min(i,1) })

  Z <- rep(0, n)
  Z[sample.int(n, n/2)] <- 1
  R <- Z * Yt + (1 - Z) * Yc

  odds.ratio <- function(y, z, b) {
    z <- as.logical(z)
    p1 <- mean(y[z])
    p2 <- mean(y[!z])

    (p1 * (1 - p2)) / (p2 * (1 - p1))
  }

  # let theta be the hypothesized value of the odds ratio
  or.moe <- function(y, z, b, theta) {
    if (theta == 1) {
      return(y)  
    }

    n <- length(y)
    nt <- sum(z) # number of treated to form column probs
    pt <- nt/n
    nr <- sum(y) # number of "correct" to form row probs
    pr <- nr/n

    # the following comes from Wikipedia's treatment of the odds ratio
    S <- sqrt((1 + (pr + pt) * (theta - 1))^2 + 4 * theta * (1 - theta) * pr * pt)

    p11 <- (1 + (pr + pt) * (theta - 1) - S)/(2 * (theta - 1))
    
    # ... 

    # or is c11 * c 22 / (c12 * c21) where c11 + c21 = nt
  }
  
  res <- parameterizedRandomizationDistribution(R, Z, test.stat = odds.ratio)

  res.fisher <- fisher.test(table(R, Z))

  res.p <- general.two.sided.p.value(res[1,1], res[1,-1])
  
  expect_equal(res.p, res.fisher$p.value)
})

