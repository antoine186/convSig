#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
NumericMatrix conv_bimap(int N, int numbase) {
  
  NumericMatrix feat_ind(numbase, N);
  int mid;
    
  if (numbase == 3) {
    mid = 1;
  } else if (numbase == 5) {
    mid = 2;
  }
  
  int temp;
    
  for (int i = 0; i < N; ++i) {
    temp = i;
    for (int j = 0; j < numbase; ++j) {
        
      if (j < mid) {
        feat_ind(j, i) = (j*4+temp%4) + 1;
        temp = floor(temp/4);
      } else if (j == mid) {
        feat_ind(j, i) = (j*4+temp%2) + 1;
        temp = floor(temp/2);
      } else if (j > mid) {
        feat_ind(j, i) = (j*4-2+temp%4) + 1;
        temp = floor(temp/4);
      }
        
    }
  }
  
  return feat_ind;
}



