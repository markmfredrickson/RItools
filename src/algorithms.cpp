#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double computeTestStatPval(Function test, NumericVector y, double t, NumericMatrix zs) {

  int sum = 0;
  int n = zs.ncol();
  
  for (int i = 0; i < n; ++i) {
    NumericVector tmp = test(y, zs(_, i));
    if (tmp[0] >= t) {
      ++sum;
    }
  }

  return (double)sum / (double)n;
}
