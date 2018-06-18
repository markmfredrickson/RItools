## Submission

This is the 0.1-16 maintenance release of RItools that addresses the problem
with undeclared package dependencies in the unit tests of the 0.1-15 release.
Our solution to this problem has been to move the offending test to a separate
file `test.notforCRAN.R` and to add that file to .Rbuildignore.

## Test environments
* local OS X install, R 3.5.0
* local Windows 10 install, R 3.4.3, x64
* win-builder (devel and release)

## R CMD check results

 * On OS X local the check was clean: 

    * R CMD check results
      0 errors | 0 warnings | 0 notes

 * On local R 3.4.3 on Windows 10 there was 1 NOTE:

   * checking package dependencies ... NOTE
     Package suggested but not available for checking: 'RSVGTipsDevice'

 * On winbuilder there was 1 NOTE:

    * checking package dependencies ... NOTE
      Package suggested but not available for checking: 'RSVGTipsDevice'

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using
`devtools::revdep_check()` and did not find any warnings, errors, or notes that
had to do with the change that we made to RItools.

