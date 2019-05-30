#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
NumericMatrix conv_bimap(int N, int K, int numbase) {
  
  NumericMatrix feat_ind(numbase, N);
  if (numbase == 3) {
    
    //int N = pow(4, (numbase-1)*2); //256
    int mid = 1;
    
    for (int i = 0; i < N; ++i) {
      int temp = i;
      for (int j = 0; j < numbase; ++j) {
        
        if (j < mid) {
          feat_ind[j, i] = j*4+temp%4;
          temp = floor(temp/4);
        } else if (j == mid) {
          feat_ind[j, i] = j*4+temp%2;
          temp = floor(temp/2);
        } else if (j > mid) {
          feat_ind[j, i] = j*4-2+temp%4;
          temp = floor(temp/4);
        }
        
      }
    }
  }
  
  return feat_ind;
}



