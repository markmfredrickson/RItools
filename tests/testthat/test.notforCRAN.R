### Tests we don't want to have run on CRAN   ###
### (ordinarily because the test will throw   ###
### a warning or error during CRAN-checking)  ###

###    tests relocated from test.utils.R      ###

test_that("data.table options issue #69", {

  if (suppressMessages(suppressWarnings(require(data.table)))) {
    data(nuclearplants)
    f <- function() 1
    expect_equal(withOptions(list(), f), 1)
  }

})
