## Test environments
* local OS X install, R 3.2.3
* win-builder (devel and release)

## R CMD check results

There was 1 NOTE with R stable from `devtools::build_win()`:

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'RSVGTipsDevice'

There were 2 NOTES with R development from `devtools::build_win()`. We are nearing release of a major revision of our package, so we are submitting this version in the expectation that we would have cleared up this NOTE before r-devel is released.

* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'RSVGTipsDevice'

* checking R code for possible problems ... NOTE
```
.balanceplot: no visible global function definition for 'axis'
balanceplot: no visible global function definition for 'par'
balanceplot: no visible global function definition for 'strwidth'
balanceplot: no visible global function definition for 'plot'
balanceplot: no visible global function definition for 'axis'
balanceplot: no visible global function definition for 'na.omit'
balanceplot: no visible global function definition for 'abline'
balanceplot: no visible global function definition for 'legend'
findStrata: no visible global function definition for 'terms'
findStrata: no visible global function definition for 'update'
findStrata: no visible global function definition for 'as.formula'
flatten.xbalresult: no visible global function definition for 'formula'
flatten.xbalresult : signifier: no visible global function definition
  for 'symnum'
makePval: no visible global function definition for 'pnorm'
naImpute: no visible binding for global variable 'median'
naImpute: no visible global function definition for 'terms.formula'
naImpute: no visible global function definition for 'get_all_vars'
naImpute: no visible global function definition for 'update.formula'
naImpute: no visible global function definition for 'as.formula'
print.xbal : <anonymous>: no visible global function definition for
  'formula'
print.xbal : <anonymous> : signifier: no visible global function
  definition for 'symnum'
print.xbal : <anonymous> : ftabler: no visible global function
  definition for 'ftable'
xBalance: no visible binding for global variable 'median'
xBalance: no visible global function definition for 'formula'
xBalance: no visible global function definition for 'terms.formula'
xBalance : <anonymous>: no visible global function definition for
  'terms'
xBalance: no visible global function definition for 'pchisq'
xBalance.find.goodstrats: no visible global function definition for
  'complete.cases'
xBalance.find.goodstrats : <anonymous> : <anonymous>: no visible global
  function definition for 'var'
xBalance.makeMM: no visible global function definition for
  'model.frame'
xBalance.makeMM: no visible binding for global variable 'na.pass'
xBalance.makeMM: no visible global function definition for 'terms'
xBalance.makepooledsd: no visible binding for global variable 'var'
xBalanceEngine : <anonymous>: no visible global function definition for
  'var'
xBalanceEngine: no visible binding for global variable 'var'
xBalanceEngine: no visible global function definition for 'pnorm'
Undefined global functions or variables:
  abline as.formula axis complete.cases formula ftable get_all_vars
  legend median model.frame na.omit na.pass par pchisq plot pnorm
  strwidth symnum terms terms.formula update update.formula var
Consider adding
  importFrom("graphics", "abline", "axis", "legend", "par", "plot",
             "strwidth")
  importFrom("stats", "as.formula", "complete.cases", "formula",
             "ftable", "get_all_vars", "median", "model.frame",
             "na.omit", "na.pass", "pchisq", "pnorm", "symnum", "terms",
             "terms.formula", "update", "update.formula", "var")
to your NAMESPACE.
```

## Downstream dependencies

I have also run R CMD check on downstream dependencies of RItools using `devtools::revdep_check()`


