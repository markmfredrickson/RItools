#include<Rcpp.h>
#include "RItools.h"

using namespace Rcpp;

// [[Rcpp::export]]
double computeTestStatPval(SEXP test, NumericVector y, double t, NumericMatrix zs) {

  int sum = 0;
  int n = zs.ncol();
  
  for (int i = 0; i < n; ++i) {
    if (callTestStatistic(test, y, zs(_, i)) >= t) {
      ++sum;
    }
  }

  return (double)sum / (double)n;
}
