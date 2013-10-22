#ifndef RITOOLS_H
#define RITOOLS_H

#include<Rcpp.h>

typedef double (*testStat)(Rcpp::NumericVector&, Rcpp::NumericVector&);
double callTestStatistic(SEXP, Rcpp::NumericVector, Rcpp::NumericVector);

#endif // RITOOLS_H

