################################################################################
# Tests for plotting functions 
################################################################################

library("testthat")

context("plot.xbal")

test_that("Issue 21: Cairo/pango errors when running plot.xbal", {
  if (capabilities()["cairo"]) {
    set.seed(20130522)

    z <- rbinom(100, size = 1, prob = 1/2)
    x <- rnorm(100)
    y <- 3 * z * x + rnorm(100)
    data <- data.frame(z, x, y)

    xb <- xBalance(z ~ x + y, data = data)
    tmpo <- tempfile()
    
    # at the moment, I haven't found a way to capture the stderr output from the C level pango function
    # so, we'll just have to know that if the errors appear in the output stream during testing
    # we should come here to see the test case (not the best strategy)

    tmpf <- tempfile()
    svg(tmpf)
    plot(xb)
    dev.off()
    file.remove(tmpf)
  }
})


