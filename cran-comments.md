## Submission

This is the 0.1-16.2 maintenance release of RItools that ensures that this
package is compatible with R 3.6.0.


## Test environments
* local OS X install, R 3.6.0
* win-builder (devel and release)

## R CMD check results

 * On OS X local the check was clean: 

    * R CMD check results
      0 errors | 0 warnings | 0 notes

 * On winbuilder there was 1 NOTE:

    * checking package dependencies ... NOTE
      Package suggested but not available for checking: 'RSVGTipsDevice'

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using
`revdepcheck::revdep_check()` and did not find any warnings, errors, or notes that
had to do with the change that we made to RItools.

