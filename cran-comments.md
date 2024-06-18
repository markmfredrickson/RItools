## Submission

This is the 0.3-4 release of RItools.  This is a maintenance release to address
errors arising from changes in `SparseM::chol()` as of SparseM-1.82

Some platforms do not like the doi formatting in DESCRIPTION although the
Checklist for R Package submissions tells us to format it exactly as we have
done: <https://cran.r-project.org/web/packages/submission_checklist.html>.

## Test environments

 - local OS X install, R 4.4.0
 - win-builder (devel and release)
 - rhub default platforms

## R CMD check results on local OS X

── R CMD check results ─────────────────────────────────────────────────────────────────── RItools 0.3-4 ────

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using
`tools::check_packages_in_dir(dir, reverse = list())`.
