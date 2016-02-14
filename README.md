Stable version: [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RItools)](http://cran.r-project.org/package=RItools)

# RItools: Randomization Inference Tools

The `RItools` package implements useful functions for implementing
randomization inference based statistical tests.  The package provides tools
for testing balance of observed covariates in observational studies using the
methodology of:

    Ben B. Hansen and Jake Bowers (2008). Covariate balance in simple,
      stratified and clustered comparative studies. Statistical Science.
      23(2):219--236.

See the online documentation for `xbalance` for more details.

The package also provides outcome analysis of simple or block randomized
trials (or matched observational studies) based on user defined models and
test statistics. See the online documentation of
`parameterizedRandomizationDistribution` for more details.

`RItools` is available on [CRAN](http://cran.r-project.org):

    > install.packages("RItools")
    > library("RItools")


##  Using a development version of RItools

These directions will install development version in a way that will not
overwrite an existing installation of `RItools` from CRAN. You will will need
to know the name of the branch you wish to install.

1. `master`: The current released version of `RItools` and a holding place for small bug changes.
2. `randomization-distribution`: Experimental work on outcome analysis using
   user defined models of effects and test statistics. This branch contains
   the tools necessary to compute estimated treatment effects, p-values, and
   confidence intervals (regions) using direct simulation from the randomization
   distribution implied by the design of the experiment (or using the exact
	   randomization distribution if the number of possible ways for the treatment to
	   be assigned is relatively small).

Install and load the `devtools` package:

    > install.packages("devtools")
    > library("devtools")

Next, pick a location to install the package. For example, create a
directory called `~/R/RItools.experimental/` (`~` is short for my home directory on a
UNIX system). For this session, we will set the library path to look in this
location first and install the package there:

    > .libPaths("~/R/RItools.experimental/") # <- your path here
    > install_github("markmfredrickson/RItools@randomization-distribution")

The function `install_github` will load the package automatically. In the
future, if you wish load the downloaded version of `RItools` in a new `R`
session you can use this one-liner:

    > library("RItools", lib.loc = "~/R/RItools.experimental") # <- your path here
