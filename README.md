[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RItools)](http://cran.r-project.org/package=RItools)

# RItools
Randomization inference tools for R

# Using a development version of RItools

## Fetching and installing in a local directory

Chances are you already have an installation of RItools that you
use. These directions will install the development version in a way
that will not overwrite your existing installation.

`ritools` is built using
[devtools](https://cran.r-project.org/package=devtools) which makes
installing the current development version very easy. Simply install
the `devtools` package and then use it to install from this
repository.

    install.packages("devtools")
    devtools:::install_github("markmfredrickson/ritools")

You may pass `ref=<branchname>` as an argument to `install_github` to
install a branch other than "master", which is the default.

Alternatively, if you have a working installation of `git` and all the
software mentioned in the previous section, you can checkout a copy of
the source directly.  Instead of downloading the source directly, fork
the project and github and clone a working copy from your forked
project:

    $ git clone git@github.com:YOURUSERNAME/ritools.git

To ensure you have all the required dependencies to work with
`ritools`, you can automatically install them with

    $ make dependencies

As mentioned, `ritools` is developed with `devtools` and requires it
to compile.  Once you have installed `devtools`, you may create a
bundled package with

    $ cd /path/to/package
    $ make build

This should build a `ritools_VERSION.tar.gz` file. You can install it
in a local directory (for example `~/R/ritools.demo`) using:

    $ mkdir -p ~/R/ritools.demo
    $ R CMD Install --no-multiarch --library=~/R/ritools.demo ./ritools_VERSION.tar.gz

You can then load the library in `R` using:

    > library("ritools", lib.loc = "~/R/ritools.demo")

### Developing for RItools

We prefer changes that include unit tests demonstrating the problem or
showing how the new feature should be added. The test suite uses the
[testthat](http://github.com/hadley/test_that) package to write and
run tests.  (Please ensure you have the latest version of testthat (or
at least v0.11.0), as older versions stored the tests in a different
directory, and may not test properly.) See the `tests/testthat`
directory for examples. To run the test suite, use:

    $ make test

New features should include inline [Roxygen](http://roxygen.org/)
documentation.  You can generate all `.Rd` documents from the
`Roxygen` code via

    $ make document

These other commands are also useful for development:

- `make interactive`: starts up an interactive session with `ritools`
  loaded.
- `make check`: runs `R CMD check` on the package
- `make vignette`: Builds any vignettes in `vignettes/` directory
- `make clean`: Removes files built by `make vignette`, `make
   document` or `make check`.  Should not be generally necessary, but
   can be useful for debugging.

When your change is ready, make a pull request on github.
