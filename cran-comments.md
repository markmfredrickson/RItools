## Resubmission

We had submitted version 0.1-14 of the package in early May, but Dr. Ligges returned 
it to us after discovering that it crashed on 64-bit Windows.  We traced the crash to 
a test of our package's interaction with the RSVGTipsDevice package that was intended
only to be run if RSVGTipsDevice was found to be installed on the system running the 
check.  This package doesn't appear to be actively maintained, and CRAN isn't 
distributing windows binaries of it; however, some systems may have incomplete 
windows installations of it, working for i386 R but not for the x64 version of R; 
perhaps CRAN's windows checking service was one of these.  In any event, our solution
is to disable the RSVGTipsDevice-invoking test for Windows platforms broadly. 

The purpose of the 0.1-14 maintenance release had been to address an issue where
our package indirectly tinkered with global options used by the data.table package.  The
0.1-15 maintenance release still addresses this issue.

## Test environments
* local OS X install, R 3.2.3
* local Windows 7 install, R 3.3.0, i386 and x64
* win-builder (devel and release)

## R CMD check results

There was 1 NOTE with R stable and development from `devtools::build_win()`:

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'RSVGTipsDevice'

## Downstream dependencies

We have also run R CMD check on downstream dependencies of RItools using `devtools::revdep_check()`


