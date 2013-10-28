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

// [[Rcpp::export]]
NumericVector runModel1(SEXP test, NumericVector y, NumericVector z, NumericMatrix zs, Function mod,  NumericVector p) {

  int np = p.length();
  double t = callTestStatistic(test, y, z);

  NumericVector result(np);

  for (int i = 0; i < np; ++i) {
    NumericVector yadj = mod(y, z, p[i]);
    result[i] = computeTestStatPval(test, yadj, t, zs);
  }

  return result;
}

// [[Rcpp::export]]
NumericVector runModel2(SEXP test, NumericVector y, NumericVector z, NumericMatrix zs, Function mod,  NumericVector p1, NumericVector p2) {

  int np1 = p1.length(), np2 = p2.length();
  double t = callTestStatistic(test, y, z);

  NumericVector result(np1 * np2);

  for (int i = 0; i < np1; ++i) {
    for (int j = 0; j < np2; ++j) {
      NumericVector yadj = mod(y, z, p1[i], p2[j]);
      result[i*np2 + j] = computeTestStatPval(test, yadj, t, zs);
    }
  }

  return result;
}

//NumericVector runModel2(Function f, NumericVector p1, NumericVector p2);
//NumericVector runModel3(Function f, NumericVector p1, NumericVector p2, NumericVector p3);
//NumericVector runModel4(Function f, NumericVector p1, NumericVector p2, NumericVector p3, NumericVector p4);
//NumericVector runModel5(Function f, NumericVector p1, NumericVector p2, NumericVector p3, NumericVector p4, NumericVector p5);
//
