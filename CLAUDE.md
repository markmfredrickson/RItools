# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

`RItools` is an R package implementing randomization-inference tools, primarily the d^2 omnibus test of Hansen and Bowers (2008, *Statist. Sci.*) for assessing covariate balance in matched observational studies and block-randomized experiments. The two user-facing entry points are `xBalance()` (older, supported but no longer the recommended interface) and `balanceTest()` (newer, supports clusters, unit weights, multiplicity adjustment, and a richer formula interface). New work should generally target `balanceTest()`; `xBalance()` is preserved for back-compatibility.

The package is on CRAN. Its current development version is tracked in `DESCRIPTION` and `NEWS.md`.

## Common commands

The `Makefile` wraps `devtools` calls in a sandboxed local library (`R_LIBS=.local`):

- `make test` --- run the full `testthat` suite via `devtools::test()`
- `make check` --- run `R CMD check` equivalent via `devtools::check()`
- `make document` --- regenerate `man/*.Rd` and `NAMESPACE` from roxygen blocks
- `make build` --- build the source tarball
- `make dependencies` --- install package dependencies
- `make interactive` --- start an R session with the package loaded via `load.R`
- `make clean` --- `git clean -Xfd` (removes all gitignored files)

Running a single test file directly (without `make`):

```sh
R -q -e 'devtools::test(filter="balanceTest")'   # filters by file name (test.balanceTest.R)
R -q -e 'devtools::test_active_file("tests/testthat/test.Design.R")'
```

The legacy `tests/*.R` scripts (e.g. `tests/balanceTestTests.R`) are run by `R CMD check` against their `.Rout.save` companions. If you intentionally change output, regenerate the `.Rout.save` files with `R CMD BATCH`.

`load.R` is a development helper used by `make interactive`: it calls `devtools::load_all(export_all = FALSE)` so that internal functions remain hidden, which matches the installed-package experience. Use this if you need to reproduce a user-visible bug locally.

GitHub Actions runs `R-CMD-check` on macOS, Windows, and Ubuntu (`devel`, `release`, `oldrel-1`) on every push and PR to `main`/`master` --- see `.github/workflows/check-standard.yaml`.

## Architecture

The package's core data flow is `formula -> Design -> engine -> xbal result -> print/plot/tidy`. Understanding the Design layer is essential before touching anything in `balanceTest.R` or `xBalance.R`.

### Design layer (`R/Design.R`)

Two S4 classes carry covariate, treatment, strata, and cluster information through the pipeline:

- `ModelMatrixPlus` --- a model matrix paired with a compact encoding of per-term missingness patterns. Its `NotMissing` slot is a numeric matrix in `[0,1]`; the first column is always all-ones (an "any X recorded" anchor), and subsequent columns describe distinct missingness patterns shared across terms. The `NM.Covariates` and `NM.terms` integer slots are look-up tables mapping covariate columns and term labels to columns of `NotMissing`. Terms with the same missingness pattern share a single column.
- `DesignOptions` (extends `ModelMatrixPlus`) --- adds `Z` (treatment), `StrataFrame` (one column per stratification, including the unstratified "`--`" column when present), and `Cluster`.

`makeDesigns(fmla, data)` is the constructor. It uses `terms(..., specials = c("cluster", "strata"))` --- relying on `survival::strata()` and `survival::cluster()`, which is why `survival` is in `Depends` (see `NEWS.md` 0.3-5 / issue #141). It validates that clusters have homogeneous treatment, that clusters nest within strata, and that each stratum contains both treated and control units.

`aggregateDesigns()` aggregates element-level designs to cluster level when clusters are present. After aggregation, `NotMissing` columns hold weighted *averages* of element-wise non-missingness, not totals --- so values can be fractional.

### Engines (`R/balanceTestEngine.R`, `R/xBalanceEngine.R`)

The engines accept the prepared design plus stratum weights and a pooled SD, then compute (per stratification):

- Adjusted within-stratum mean differences combined across strata (ETT-type weighting for descriptives; harmonic-times-mean-weight by default for inferentials).
- Univariate z-scores via the permutation null variance, with `p.adjust()` applied (`balanceTest` defaults to Holm).
- The omnibus chi-square via SVD of the (treatment-centered, weight-scaled) covariate matrix --- the d^2 statistic of Hansen and Bowers (2008).
- A `tcov` test-statistic covariance matrix.

`post.alignment.transform` (e.g. `rank` for a Wilcoxon-like test) operates on the stratum-aligned `tmat` *after* centering and *before* the variance/SVD calculations; the engine recenters after the transform.

### Stratum weighting

`harmonic_times_mean_weight()` (in `R/harmonic.R`) is the default `stratum.weights` for `balanceTest()` inferentials. It multiplies the per-stratum harmonic mean of treated/control counts by the per-stratum mean of `unit.weights`, which is optimal under the modeling assumptions of Kalton (1968) and Hansen and Bowers (2008, secs. 3.2 and 5). Descriptive (mean-difference) calculations in `balanceTest()` use a different, ETT-type weighting; do not conflate the two.

### Result objects

`balanceTest()` returns an object of class `c("balancetest", "xbal", "list")`; `xBalance()` returns `"xbal"`. Most methods (`print`, `plot`, `xtable`, `subset`, `tidy`, `glance`) live on `"xbal"` and are inherited. The 3-D `results` array is indexed `[variable, statistic, stratification]`; the unstratified comparison uses the stratification name `"--"`.

### NA handling

Missingness is handled in two different ways depending on which calculation is being performed:

- **Descriptives** drop missing values.
- **Inferentials** impute to within-stratum means (see `naImpute.R`) so that the test statistic and its null variance are computed on a complete matrix.

Item-level missingness proportions are reported in the printed output, parenthesized and pushed to the bottom of the variable list. Variables that share a missingness pattern share a single missingness column, so the displayed missingness label may not match the variable name --- this is intentional, not a bug.

## Conventions worth knowing

- Documentation is roxygen2 (see `RoxygenNote` in `DESCRIPTION`); `man/*.Rd` and `NAMESPACE` are generated --- edit the roxygen blocks, then `make document`.
- `lexicon.txt` + `checkspelling.R` provide an aspell-based spellcheck for `.Rd` files. Add genuinely new technical terms to `lexicon.txt` rather than rewording prose.
- `tests/testthat/test.notforCRAN.R` is excluded from CRAN builds via `.Rbuildignore`. Tests that require optional packages, are slow, or rely on external resources should go there, not in the main test files.
- The legacy `tests/*.R` + `*.Rout.save` pairs predate `testthat` and are still used by `R CMD check`. Prefer adding new tests under `tests/testthat/`.
- Internal helpers are kept unexported deliberately; `load.R` mirrors this with `export_all = FALSE`. Do not add `@export` to internal helpers without a clear reason.
- `cran-comments.md` and `CRAN-SUBMISSION` are part of the CRAN release process; the `revdep/` directory holds reverse-dependency check output. Leave these alone unless preparing a release.
