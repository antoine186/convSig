#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
RcppExport SEXP RM_nonSNP(SEXP x) {
  NumericMatrix startend = as<NumericMatrix>(x);
    
  int nr = startend.nrow();
  int nc = startend.ncol();
  
  int* rm_ind = new int[nr];
  
  for (int i = 0; i < nr; ++i)  
  { 
    if (startend(i,1) == startend(i,2)) {
      rm_ind[i] = 1;
    } else {
      rm_ind[i] = 0;
    }
  } 
  
  return wrap(rm_ind);
}

