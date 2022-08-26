approxEq <- function(a, b, e=1e-10) {
    abs(a - b) < e
}


Rcpp::sourceCpp(code = '
#include "Rcpp.h"

// [[Rcpp::export]]
double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}
')


Rcpp::sourceCpp(code = '
#include "Rcpp.h"

// [[Rcpp::export]]
double l2(Rcpp::NumericVector a) {
    return sqrt(sum(pow(a, 2.0)));
}
')


Rcpp::sourceCpp(code = '
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double dot(arma::vec x, arma::vec y) {
    return arma::dot(x, y);
}
')


Rcpp::sourceCpp(code = '
#include "Rcpp.h"

// [[Rcpp::export]]
double lgammafn(double x) {
    return R::lgammafn(x);
}
')


