################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("Plotting Functions")

test_that("Basic plot", {
  # create a quick xBalance result
  set.seed(20121119)
  Z <- rep(c(0,1), 10)
  df <- data.frame(Z = Z,
                   X1 = rnorm(20, mean = Z*2),
                   X2 = rnorm(20, mean = Z*3),
                   X3 = rnorm(20, mean = Z * -1),
                   X4 = rnorm(20, mean = Z * -0.5))

  xb <- xBalance(Z ~ ., data = df)

  x11() # this might be more common than quartz
  expect_true(dev.capabilities()$capture)

  plot(xb)
  p1 <- dev.capture()
  dev.off()

  x11()
  plot(xb)
  p2 <- dev.capture()
  expect_true(class(p1) == "matrix" & class(p2) == "matrix")
  expect_identical(p1, p2) # just to prove that the same plot twice is really identical
  dev.off()

  opts <- options(warn = 2)
  # has a an argument to make only a right-sided abs difference plot
  # it will be a warning if absolute isn't a parameter
  x11()
  plot(xb, absolute = T)

  p2 <- dev.capture()
  expect_true(!identical(p1,p2))
  dev.off()

  options(opts)

  x11()
  # has an argument to order the variables (from bottom to bottom)
  plot(subset(xb, vars = c("X1", "X2", "X4", "X3")), ordered = F)
  expect_true(!identical(p1, dev.capture()))
  dev.off()
 
  # the order the data based on the selected variable
  x11()
  plot(xb, ordered = T)
  expect_true(!identical(p1, dev.capture()))
  dev.off()
  # note: the order should change when absolute = T, but we can't really test it as the plot is then different for two reasons
  # one way to test this would be to have a helper function that creates an array for something else to plot -- and the intermediate data could be checked

  # just a sanity check to make sure that the previous dev.capture tests worked
  x11()
  plot(xb)
  expect_identical(p1, dev.capture())
  dev.off()

  ### Error checking
  expect_error(plot(xb, statistic = "foo"), "statistic")
  expect_error(plot(xb, variable.labels = c("foo")), "labels")
  expect_error(plot(xb, strata.labels = c("foo")), "labels")
})

test_that("Generic balance plots", {
  # the balanceplot function takes a matrix and plots it in the expected fashion
  # this serves as a compliment to the .xbal.plot function

  testmat <- matrix(c(4,3,2,1,-5, 3,-2,-3,2,3), ncol = 2, 
                    dimnames = list(c("Variable 1","Variable Two","Var3","X4", "V5"),
                                    c("Stratification 1", "Stratification 2")))
  grps <- as.factor(c("Grp1", "Grp2", "Grp2", "Grp1", "Grp1"))

  x11() # this might be more common than quartz
  expect_true(dev.capabilities()$capture)

  balanceplot(testmat)
  p1 <- dev.capture()
  dev.off()
  
  x11()
  balanceplot(testmat) 
  p2 <- dev.capture()
  dev.off()

  expect_identical(p1,p2)

  x11()
  balanceplot(testmat, grps)
  p.grps <- dev.capture()
  dev.off()

  expect_true(!identical(p1, p.grps))

  grps2 <- grps
  grps2[1] <- "Grp2"

  x11()
  balanceplot(testmat, grps2)
  expect_true(!identical(p.grps, dev.capture()))
  dev.off()
})
