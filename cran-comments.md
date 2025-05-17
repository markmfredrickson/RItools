## Submission

This is the 0.3-5 release of RItools.  This is a maintenance release to address
errors arising from changes in r-devel having to do with problems with
stats:::terms.formula and 'specials' and the survival package.

Some platforms do not like the doi formatting in DESCRIPTION although the
Checklist for R Package submissions tells us to format it exactly as we have
done: <https://cran.r-project.org/web/packages/submission_checklist.html>.

## Test environments

- local OS X install, R 4.5.0
- win-builder (devel and release)
- rhub default platforms

## R CMD check results on local OS X

── R CMD check results ───────────────────────────────────────────────────────────────────────────────────────────────────────────── RItools 0.3-5 ────
Duration: 38.6s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking dependencies in R code ... NOTE
  Package in Depends field not imported from: ‘survival’
    These packages need to be imported from (in the NAMESPACE file)
    for when this namespace is loaded but not attached.

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using
`tools::check_packages_in_dir(dir, reverse = list())`.
