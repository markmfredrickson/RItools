Stable version: [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RItools)](https://CRAN.R-project.org/package=RItools)

[![Travis-CI Build Status](https://travis-ci.org/markmfredrickson/ritools.svg?branch=master)](https://travis-ci.org/markmfredrickson/ritools)

# RItools: Randomization Inference Tools

The `RItools` package implements useful functions for implementing
randomization inference based statistical tests.  The package provides tools
for testing balance of observed covariates in observational studies using the
methodology of:

    Ben B. Hansen and Jake Bowers (2008). Covariate balance in simple,
      stratified and clustered comparative studies. Statistical Science.
      23(2):219--236.

See the online documentation for `xBalance` for more details.

The package also provides outcome analysis of simple or block randomized
trials (or matched observational studies) based on user defined models and
test statistics. See the online documentation of
`parameterizedRandomizationDistribution` for more details.

`RItools` is available on [CRAN](https://CRAN.R-project.org):

    > install.packages("RItools")
    > library("RItools")


##  Using a development version of RItools

These directions will install development version in a way that will not
overwrite an existing installation of `RItools` from CRAN. You will will need to
know the name of the branch you wish to install.

1. `master`: The current released version of `RItools` and a holding place for
   small bug changes.
2. `randomization-distribution`: Experimental work on outcome analysis using
   user defined models of effects and test statistics. This branch contains the
   tools necessary to compute estimated treatment effects, p-values, and
   confidence intervals (regions) using direct simulation from the randomization
   distribution implied by the design of the experiment (or using the exact
   randomization distribution if the number of possible ways for the treatment
   to be assigned is relatively small).

We recommend using `dev_mode` from the `devtools` package to install
in-development versions of the package so that you can keep the current CRAN
version as the primary package.
Activating `dev_mode` creates a secondary library of packages which can
only be accessed while in `dev_mode`. Packages normally installed can still be
used, but if different versions are installed normally and in `dev_mode`, the
`dev_mode` version takes precedent if in `dev_mode`.

Install and load the `devtools` package:

    > install.packages("devtools")
    > library("devtools")

Activate `dev_mode`:

    > dev_mode()
    d>

Note that the prompt changes from `>` to `d>` to let you know you're in
`dev_mode`. Now choose the development branch you want to use. To install
`master`:

    d> install_github("markmfredrickson/RItools")

or to install the `randomization-distribution` branch:

    d> install_github("markmfredrickson/RItools@randomization-distribution")

Either way, the package is then loaded in the usual fashion, provided
you're still in `dev_mode`:

    d> library(RItools)

Once you've done this you can disable `dev_mode` as follows

    d> dev_mode()
    >

The development version of the package remains loaded.

Note that if you load the package -- ie, enter `library(RItools)` (when the package
hasn't already been loaded otherwise) -- while _not_ in `dev_mode`,
then you'll get whatever version of the package may be installed in
your library tree, not this development version.

If you want to switch between versions of `RItools`, we suggest re-starting R.
