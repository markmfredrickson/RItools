#include<Rcpp.h>
#include "RItools.h"

using namespace Rcpp;

// [[Rcpp::export]]
double wilcoxTestStatistic(NumericVector& y, NumericVector& z) {
  int nx = sum(z);
  int ny = z.length() - nx;
  
  // it would be nice to get direct access to do_rank from R's internals.
  Environment base("package:base");
  Function rank = base["rank"];
  NumericVector r = rank(y);

  double w = sum(ifelse(z == 1, r, 0)) - (nx * (nx + 1)/2);
  
  return w;
}

// [[Rcpp::export]]
double meanDifference(NumericVector& y, NumericVector& z) {

  int len = z.length();
  int nt = sum(z);
  int nc = z.length() - nt;
  double sumt = 0.0;
  double sumc = 0.0;
 

  for(int i = 0; i < len; ++i) {
    if (z[i] == 1) {
      sumt += y[i]; 
    } else {
      sumc += y[i];
    }
  }

  return sumt/nt - sumc/nc;
}

// [[Rcpp::export]]
XPtr<testStat> testStatisticPtr(std::string f) {

  if (f == "wilcoxTestStatistic") {
    return(XPtr<testStat>(new testStat(&wilcoxTestStatistic)));
  } else if (f == "meanDifference") {
    return(XPtr<testStat>(new testStat(&meanDifference)));
  }

  return XPtr<testStat>(R_NilValue);
}

// [[Rcpp::export]]
double callTestStatistic(SEXP f_, NumericVector y, NumericVector z) {
  if (TYPEOF(f_) == EXTPTRSXP) {
    XPtr<testStat> tsp(f_);
    testStat f = *tsp;
    return f(y, z);
  }
  // should probably test this more explicitly
  Function f(f_);
  NumericVector a = f(y, z);
  return a[0];
}
