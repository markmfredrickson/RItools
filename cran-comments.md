## Submission

This is the 0.1-18 maintenance release of RItools that ensures that this
package is compatible with R 4.1.2.


## Test environments
* local OS X install, R 4.1.2
* win-builder (devel and release)

## R CMD check results

 * On OS X local there was 1 NOTE:

    * R CMD check results
      0 errors | 0 warnings | 1 note
      Package suggested but not available for checking: ‘RSVGTipsDevice’

 * On winbuilder there was 1 NOTE:

    * checking package dependencies ... NOTE
      Package suggested but not available for checking: 'RSVGTipsDevice'

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using
`revdepcheck::revdep_check()`:

### revdepcheck results

We checked 3 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
