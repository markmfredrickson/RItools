# RItools {NEXTVERSION}

* Bug fix (base R plot method svg export)
* More graceful handling in balanceTest() of the edge case of balance testing for all-constant variables, with messages not errors and p-values of NA not 1

# RItools 0.3-5.9004 (devel branch: sigma_x omnibus)

This is a development version on the `devel-sigma-x-omnibus` branch.  It is
not intended to replace the released version on `main`.

* This version merges `main` (pull request #143).  When a stratification has no covariate
  left to test (every covariate constant overall or within strata, so the
  omnibus has df 0), `balanceTest()` prints one message, containing main's
  phrase "perhaps all covariates are constant (within strata)?", instead of
  the branch's two earlier messages.
* The Cauchy combination (`cauchy.combination = TRUE`) now uses the truncated
  conversion of Gui, Jiang, and Wang (2025, Biometrika 112(4), asaf038), which
  the riposte package adopted in 0.0.0.9008.  Each per-covariate p-value p
  becomes tan((0.5 - 0.9 p) pi).  With n p-values whose converted values sum
  to S, `cauchy_comb_p` is min(1, n P(C > S) / 0.9) for a standard Cauchy
  variable C.  Liu and Xie's conversion, tan((0.5 - p) pi), sends p = 1 to
  minus infinity.  So when one covariate's adjusted difference was exactly 0,
  `cauchy_comb_p` was about 1: 0.99999999999999933 in a 40-unit example where
  the truncated version gives 4.6e-05.  When the null hypothesis is true, the
  truncated version gives a p-value below 0.05 more than 5 percent of the
  time (5.7 percent with 6 independent covariates, 6.3 percent with 40).
  Below 0.01 it stays near 1 percent.  See
  `vignettes/truncated-cauchy-memo.md`.

* New optional alternative omnibus test in `balanceTest()`: when called with
  `sigma_x_test = TRUE`, the result table `$overall` is augmented with three
  columns (`sigma_x`, `sigma_x.df`, `sigma_x.p.value`) reporting
  `T = d' Sigma_x^{-1} d` --- a Mahalanobis-type statistic that standardizes
  the adjusted differences vector by the within-stratum-pooled sample
  covariance of the covariates rather than by the permutation covariance
  used by the existing Hansen-Bowers d^2 test.  Both stats are reported
  side-by-side; the d^2 columns are unchanged.
* Four null-distribution backends are exposed via the new `null` argument
  to `balanceTest()`:
    * `"satterthwaite_finite"` (default): Satterthwaite a*chi^2_v match using
      the *exact* moments of T under the within-stratum randomization
      distribution.  Honest at finite n.  The second-order moment machinery
      is ported from work by Mark Fredrickson on the `i113-highermoments`
      branch.
    * `"satterthwaite_asymptotic"`: same Satterthwaite match using the
      Gaussian identities `sum lambda_k` and `2 sum lambda_k^2`.  Conservative
      at small n; collapses to the finite-sample backend in the large-n limit.
    * `"imhof"` and `"davies"`: analytic CDF inversion via the `CompQuadForm`
      package (now in Suggests:).
    * `"simulate"`: direct Monte Carlo over the within-stratum randomization
      distribution; configurable via the `n_simulate` argument.
* New helper `default_sigma_x()` (internal) computes the default
  within-stratum-pooled sample covariance.  Users may also supply their own
  `sigma_x` matrix via the new argument; the test stat and its null
  distribution are invariant under positive rescaling of `sigma_x`.
* New optional omnibus in `balanceTest()`: with `cauchy.combination = TRUE`,
  `$overall` gains a `cauchy_comb_p` column reporting the Cauchy combination
  (ACAT; Liu and Xie 2020) of the per-covariate p-values.  Because it combines
  the marginal tests instead of inverting their covariance, it stays
  informative when the number of covariates reaches or exceeds the rank of the
  d^2 covariance --- the high-dimensional regime in which the chi-square
  omnibus degenerates (its statistic freezes at its rank and its p-value parks
  near 0.46).  It combines the raw per-covariate p-values, missingness-
  indicator comparisons included.  Off by default; the default `$overall`
  table is unchanged.

# RItools 0.3-5

* Now "Depends" on the `survival` package to bring in the `strata` and `cluster`
functions (#141).

# RItools 0.3-4

* Using `.lm.fit()` within `slm.fit.csr.fixed` given breaking changes in `SparseM::chol()` (#134)

# RItools 0.3-3

* Adjusted arguments to `xtable.xbal()` to match the generic. (#131)

# RItools 0.3-2

* Adjusted terminology on `balanceplot()` plots when using `absolute=TRUE`.
  (#129)
* Updated some **tidyverse**/**ggplot2** code to deal with deprecations.
* Using `bibentry()` in citation to address deprecation.


# RItools 0.3-1

* Skipped tests that cause problems with systems using the ATLAS linear algebra
  library.

# RItools 0.3-0

* Exporting `broom::glance()` and `broom::tidy()` methods for balance test
  objects, a work in progress (#90).
* Addressed error relating to setting of pseudo-inversion tolerances (#103).
  xBalance stops short of correcting the error, to maintain
  back-compatibility; but it gains a new argument, pseudoinversion_tol, with
  which it can readily be fixed.
* New balanceTest function extends functionality of xBalance:
  * Support for clustered designs
  * More use of formula interface
  * Multiple comparison adjustments for individual covariate tests
  * Support for unit level weighting
* Introduced ggplot2 based plotting. Optional for xBalance objects; default for
  new balanceTest objects.

# RItools 0.2

This was an internal release not placed on CRAN.

# RItools 0.1-18

This is a maintenance release.

* Ensured that the package was compatible with R 4.1.2

# RItools 0.1-17

This is a maintenance release.

* Ensured that the package was compatible with R 3.6.0

# RItools 0.1-16

This is a maintenance release.

* Fixed a problem where the test suite was calling an external package that was
  not a formal dependency of the core functionality of the package.

# RItools 0.1-15

This is a maintenance release.

* Fixed an obscure bug in which we indirectly fiddled with `data.table` related
  global options (#69)
* Disabled tests of RSVGTipsDevice-dependent functionality on Windows platforms,
  where (as of this writing) RSVGTipsDevice does not build and check reliably
  and is not distributed in binary through CRAN (#71)

# RItools 0.1-14

This version was submitted to CRAN prematurely and ultimately was not released.

# RItools 0.1-13

* With `xBalance()`, you can now specify strata "foo" and "bar" by including "+
  strata(foo) + strata(bar)" in the `fmla` argument, without need to give a
  separate `strata` argument (eg `strata=list(foo=~foo, bar=~bar)`).

# RItools 0.1-12

* Switched to base graphics instead of the lattice package for balance
  plots.
* Option to make plots on absolute, instead of signed, values.
* Added the ability to include tooltips to balance plots when using
  the RSVGTipsDevice package.
* Added a `balanceplot()` function to handle plotting matrices, not just `xbal`
  objects.
* Balance plots can be grouped into related variables. By default
  factors are automatically grouped for `xbal` objects.
* Some performance improvements when certain results are not requested.
* `post.alignment.transform` argument to `xBalance()` allows modifying the data
  after per-stratum centering.
* `subset` method for `xbal` objects.

# CHANGES in RItools v. 0.1-11

## NEW FEATURES
* the null SD of the difference in adjusted means is now reportable as
  "adj.means.diffs.sd"
* `report="all"` provides all relevant statistics
* Entering `args(xBalance)` is now useful as a reminder of report= options.

## DEPRECATED & DEFUNCT
* The default imputation method for missing data is the median. Before 0.1-9 it
  was the mean. To use the mean, use `impfn=mean.default`.
* Passing `strata=NULL` to `xBalance()` is depracated, as it can be easy to do
  inadvertently, with a misspelling of a non-null intended argument. For no
  stratification, use `list(unstrat=NULL)`, as is now the default.

##  BUG FIXES
* `display`, `align`, etc arguments to `xtable.xbal()` have been enabled.
* In `xBalance()`, when specifying a stratification using a formula with several
  variables, as in `~x+y`, the stratification you get is now interaction(x,y)
  rather than just y.
