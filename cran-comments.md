## Submission

This is the 0.3-2 release of RItools.  This is a maintenance release.
We updated some **tidyverse**/**ggplot2** code to deal with deprecations.

Some platforms do not like the doi formatting in DESCRIPTION although the
Checklist for R Package submissions tells us to format it exactly as we have
done: <https://cran.r-project.org/web/packages/submission_checklist.html>.

## Test environments

 - local OS X install, R 4.2.2
 - win-builder (devel and release)
 - rhub default platforms

## R CMD check results on local OS X

── R CMD check results ─────────────────────────────────────────────────────────────────── RItools 0.3-2 ────

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using
`revdepcheck::revdep_check()`:
