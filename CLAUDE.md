# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project overview

RItools (Randomization Inference Tools) is an R package for covariate balance
assessment in simple, stratified and clustered comparative studies, after
Hansen and Bowers (2008, Statistical Science). The current entry point is
`balanceTest()`; `xBalance()` is the older interface it replaces and is kept
for back-compatibility.

This clone is the upstream, `markmfredrickson/RItools` (default branch
`main`). MEERKT (the `ee_statistical_component` repository) installs RItools
from the fork `benbhansen-stats/RItools`, which is read only for us. Work here
with the understanding that changes are ultimately meant for MEERKT and reach
it only when that fork is updated.

## Build and test

Local R lacks most of the package's dependencies. The MEERKT image
`meerkt/rpkgs` has them all (devtools, testthat, roxygen2, SparseM, svd,
optmatch, ...), so run everything in that image with the clone mounted:

```bash
# whole suite (about a minute; 1506 tests pass on main at 794ebcc)
docker run --rm -v "$PWD":/pkg -w /pkg meerkt/rpkgs:latest \
  Rscript -e 'devtools::test()'

# one file
docker run --rm -v "$PWD":/pkg -w /pkg meerkt/rpkgs:latest \
  Rscript -e 'devtools::test(filter = "balanceTest")'

# regenerate NAMESPACE and man/ after changing roxygen comments
docker run --rm -v "$PWD":/pkg -w /pkg meerkt/rpkgs:latest \
  Rscript -e 'devtools::document()'

# R CMD check
docker run --rm -v "$PWD":/pkg -w /pkg meerkt/rpkgs:latest \
  Rscript -e 'devtools::check()'
```

Use `devtools::test()`, not `load_all()` plus `test_dir()`: the latter skips
S4 method registration and reports spurious failures in `test.Design.R` and
`test.tidiers.R`.

A quirk of the formula method, present before the generic: `unit.weights`
and `subset` go through `model.frame()` without the formula, so an
expression there is evaluated against the data frame and then the global
environment, not the caller's frame. Inside a function or a test, pass a
column of the data frame, not an expression built from local objects.

The `Makefile` targets (`make test`, `make check`, `make document`) wrap the
same devtools calls but assume a local `.local` library; they do not work on
this machine without the image.

Tests live in `tests/testthat/`. The top-level `tests/*.R` with `.Rout.save`
files are older regression scripts run by `R CMD check`. `.Rbuildignore`
keeps `test.notforCRAN.R` and `test.tidiers.R` out of the CRAN build.

## Architecture

`balanceTest()` is an S3 generic (`R/balanceTest.R`) with two methods that
share one internal body, `.balanceTestFromFrame()`, which runs steps 2 to 8
below from a formula and a model frame carrying a `(weights)` column:

- `balanceTest.formula`: the original interface, with its first formal still
  named `fmla` so calls by name keep working. Step 1 below is its own. Both
  methods stop on an argument they do not recognize (`.stopOnUnusedArguments`).
- `balanceTest.balanceSpec` (`R/balanceSpec.R`): a `balanceSpec()` object
  carries treatment, strata, optional cluster and unit weights, the variables
  to test, and a supplied covariance matrix of the expanded covariates, whose
  pseudo-inverse square root is computed once in the constructor. The method
  builds the formula, sets `covariate.scales` from the positive entries of
  the matrix's diagonal, and runs the shared body with `HB08_fixed(spec)` as
  the inferentials calculator: the null covariance of the adjusted
  differences is the supplied matrix times a design scalar (sum over strata
  of the squared stratum weight ratio over `1/n_t + 1/n_c`), so the inverse
  is reused. Where values were imputed, each entry is instead scaled stratum
  by stratum by the share of units observed on both variables (read from the
  new `NotMissing` slot of `CovsAlignedToADesign`), and the inverse is
  recomputed. Expanded covariate columns missing from the matrix are an
  error; matrix columns the data do not produce are set aside with one
  message and the inverse for the remaining block is cached in the closure.
- `expandCovariates(fmla, data)` (exported, `R/Design.R`) returns the
  expanded covariate columns exactly as `makeDesigns()` makes them (one
  indicator per factor level, character columns as factors, NAs kept), via
  the shared internal `.expandModelFrame()`. It is what a caller computing a
  covariance elsewhere should use so the columns line up.

The `balanceTest()` pipeline:

1. `model.frame` on the formula with `unit.weights` mapped to `(weights)` and
   `subset` smuggled through the `offset` slot (hence no `...` in the formula).
2. `makeDesigns(fmla, data)` (`R/Design.R`): parses `strata()` and
   `cluster()` specials from the right-hand side, checks that clusters are
   nested in strata and carry one treatment value, drops strata without both
   groups (warning), expands factors with one dummy per level, and returns a
   `DesignOptions` object (S4; extends `ModelMatrixPlus`, which carries the
   covariate matrix plus a `NotMissing` matrix of missingness patterns).
   The unstratified comparison is the stratification named `--`, present
   unless the formula has `-1`.
3. `aggregateDesigns()` collapses elements to clusters.
4. `designToDescriptives()` computes the descriptive columns (`Control`,
   `Treatment`, `adj.diff`, `std.diff`, `pooled.sd`) with effect of treatment
   on treated stratum weights; `stratum.weights` does not affect these.
5. `DesignWeights(aggDesign, stratum.weights)` produces per-stratum
   `sweights` (normalized to sum to 1) and `wtratio` for the inferential
   calculations. Default `harmonic_times_mean_weight` is the harmonic mean of
   treated and control counts times the stratum mean of unit weights. Note
   the factor of 2 relative to the paper's h_b, documented inline.
6. `alignDesignsByStrata()` centers covariates within strata, imputes missing
   values to the stratum mean, applies `post.alignment.transform`, and
   returns a `CovsAlignedToADesign`.
7. `HB08()` (the `inferentials.calculator`) computes the adjusted differences
   `ssn`, their permutation covariance `tcov`, univariate z and p, and the
   omnibus chi-square `Msq` with `DF` from `XtX_pseudoinv_sqrt`
   (`R/utils.R`). Zero-variance covariates get z and p of NA and are left out
   of the omnibus statistic. `check_for_degenerate()` warns when the rank is
   at least the number of units minus strata. `HB08_2016` is the earlier
   version, kept for comparison.
8. Back in `balanceTest()`: p-values are adjusted per stratification with
   `p.adjust(method = p.adjust.method)` (default `holm`); not-missing
   indicators that are constant are removed; the result is a list of class
   `c("balancetest", "xbal", "list")` with `$overall` (chisquare, df, p.value
   per stratification, with `tcov` as an attribute) and `$results` (a
   vars x stat x strata array with attributes `NMpatterns`, `originals`,
   `term.labels`).

`xBalance()` and `balanceTestEngine()` are the old path; `print.xbal`,
`plot.balancetest`, `xtable.xbal`, `flatten.xbalresult` and `xbal_tidiers`
(broom `tidy`/`glance`) consume the result object.

## How MEERKT uses this package

Mostly `balanceTest()` on a matched design (`z ~ covs + strata(match)`), and
`RItools:::flatten.xbalresult()` for HTML output. `permutationAnalysis()` in
meerktjob `analysis.R` also uses it as a permutation test engine for outcomes:
`z ~ y + strata(match) - 1` with a `post.alignment.transform` that trims, and
reads `$results` and the `tcov` attribute of `$overall`. Changes to the shape
of the result object or to the treatment of zero-variance and degenerate cases
reach MEERKT directly.

## Conventions

- roxygen comments generate `NAMESPACE` and `man/`; run `devtools::document()`
  after editing them and commit the output.
- `NEWS.md` has a `{NEXTVERSION}` section for unreleased changes; add a line
  there for user-visible changes.
- Commit subjects in this repo often carry a prefix (`BUG`, `ENH`, `DOC`,
  `TST`, `STY`).
- Do not name an internal function `balanceTest.<something>`: with
  `balanceTest` a generic, R CMD check reads it as an unregistered method.
- Open issues of note: #142 (Cauchy combination to side-step degeneracy, on
  branch `cauchycomb`), #25 (flagging degenerate chi-square distributions),
  #113 (higher moments, branch `i113-highermoments`), #90 (broom tidiers),
  #34 (p.adjust for univariate p-values; now implemented).
