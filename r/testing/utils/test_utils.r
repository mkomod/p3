library(Rcpp)

Rcpp::sourceCpp("../src/utils.cpp")
s = rep(1:5, each=3)

get_group_indices(s)
