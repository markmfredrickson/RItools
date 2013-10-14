#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double wilcoxTestStatistic(NumericVector y, NumericVector z) {
  int nx = sum(z);
  int ny = z.length() - nx;
  
  // it would be nice to get direct access to do_rank from R's internals.
  Environment base("package:base");
  Function rank = base["rank"];
  NumericVector r = rank(y);

  double w = sum(ifelse(z == 1, r, 0)) - (nx * (nx + 1)/2);
  
  return w;
}
