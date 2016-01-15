## Test environments
* local OS X install, R 3.2.3
* win-builder (devel and release)

## R CMD check results

There was 1 NOTE with R stable and development from `devtools::build_win()`:

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'RSVGTipsDevice'

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using `devtools::revdep_check()`


