#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// This is a for loop comparison helping to remove non single nucleotide changes

// [[Rcpp::export]]
LogicalVector RM_nonSNP(DataFrame startend, SEXP ar) {
  LogicalVector rm_ind = as<LogicalVector>(ar);
  
  NumericVector startend1 = startend[0];
  NumericVector startend2 = startend[1];
  
  int nr = startend.nrow();
  
  for (int i = 0; i < nr; ++i)  
  { 
  if (startend1[i] == startend2[i]) {
  rm_ind[i] = TRUE;
  } else {
  rm_ind[i] = FALSE;
  }
  } 
  
  return rm_ind;
}
