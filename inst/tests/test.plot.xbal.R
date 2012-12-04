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
  plot(xb)
  p2 <- dev.capture()
  expect_true(class(p1) == "matrix" & class(p2) == "matrix")
  expect_identical(p1, p2) # just to prove that the same plot twice is really identical

  opts <- options(warn = 2)
  # has a an argument to make only a right-sided abs difference plot
  # it will be a warning if absolute isn't a parameter
  plot(xb, absolute = T)

  p2 <- dev.capture()
  expect_true(!identical(p1,p2))

  options(opts)

  # has an argument to order the variables (from bottom to bottom)
  plot(xb, which.vars = c("X1", "X2", "X4", "X3"))
  expect_true(!identical(p1, dev.capture()))
 
  # but including other variables is an error
  expect_error(plot(xb, which.vars = c("X1", "X2", "X3", "X4", "X5")), "Unknown variable\\(s\\): X5")

  # the order the data based on the selected variable
  plot(xb, ordered = T)
  expect_true(!identical(p1, dev.capture()))
  # note: the order should change when absolute = T, but we can't really test it as the plot is then different for two reasons
  # one way to test this would be to have a helper function that creates an array for something else to plot -- and the intermediate data could be checked

  # just a sanity check to make sure that the previous dev.capture tests worked
  plot(xb)
  expect_identical(p1, dev.capture())

  dev.off()
})

test_that("Helper function", {
  set.seed(20121119)
  Z <- rep(c(0,1), 10)
  df <- data.frame(Z = Z,
                   X1 = rnorm(20, mean = Z*2),
                   X2 = rnorm(20, mean = Z*3),
                   X3 = rnorm(20, mean = Z * -1),
                   X4 = rnorm(20, mean = Z * -0.5))

  xb <- xBalance(Z ~ ., data = df, strata = data.frame(raw = factor('none'), other = as.factor(rep(c(1,0), each = 10))))

  lbls <- c("X One", "X Two", "X Three", "X Four")
  res <-  .plot.xbal(xb, 
                    which.strata = c("raw", "other"),
                    thestratalabs = c("Unstratified", 
                                      "Some Other Strata"),
                    which.stat = "std.diff",
                    which.vars = c("X1", "X2", "X3", "X4"),
                    thevarlabs = lbls ,
                    absolute = F,
                    ordered = F)

  expect_equal(dim(res), c(4,2))
  expect_equal(rownames(res), lbls)
  expect_equal(colnames(res), c("Unstratified", "Some Other Strata"))

  res.ordered <-  .plot.xbal(xb, 
                    which.strata = c("raw", "other"),
                    thestratalabs = c("Unstratified", 
                                      "Some Other Strata"),
                    which.stat = "std.diff",
                    which.vars = c("X1", "X2", "X3", "X4"),
                    thevarlabs = lbls ,
                    absolute = F,
                    ordered = T)

  expect_false(identical(res, res.ordered))
  expect_equal(dim(res.ordered), c(4,2))
  expect_equal(rownames(res.ordered), lbls[c(3,4,1,2)])

  res.pos <-  .plot.xbal(xb, 
                    which.strata = c("raw", "other"),
                    thestratalabs = c("Unstratified", 
                                      "Some Other Strata"),
                    which.stat = "std.diff",
                    which.vars = c("X1", "X2", "X3", "X4"),
                    thevarlabs = lbls ,
                    absolute = T,
                    ordered = T)

  expect_false(identical(res.pos, res.ordered))
  expect_equal(dim(res.pos), c(4,2))
  # the row order of the final two variables should flip when taking abs
  expect_equal(rownames(res.pos), lbls[c(4,3,1,2)])
  expect_true(all(res.pos >= 0))

  # repeating with a single strata
  
  res.singlestrat <-  .plot.xbal(xb, 
                    which.strata = c("raw"),
                    thestratalabs = c("Unstratified"),
                    which.stat = "std.diff",
                    which.vars = c("X1", "X2", "X3", "X4"),
                    thevarlabs = lbls ,
                    absolute = F,
                    ordered = T)

  expect_equal(dim(res.singlestrat), c(4,1)) # should be a column data.frame, not a vector

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
  balanceplot(testmat) 
  p2 <- dev.capture()
  expect_identical(p1,p2)

  balanceplot(testmat, grps)
  p.grps <- dev.capture()
  expect_true(!identical(p1, p.grps))

  grps2 <- grps
  grps2[1] <- "Grp2"
  balanceplot(testmat, grps2)
  expect_true(!identical(p.grps, dev.capture()))

  dev.off()
})
