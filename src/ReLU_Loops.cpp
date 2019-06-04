#include <Rcpp.h>
#include <vector>
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

// // [[Rcpp::export]]
// int four_rec(NumericVector flat_ar, NumericVector g_dims, int ax, int want_len) {
//   
//   int prelim_zero = g_dims[0];
//   int prelim_one = g_dims[1];
//   int prelim_two = g_dims[2];
//   int prelim_three = g_dims[3];
//   
//   // Reconstruct flat_ar
//   //double* recon_ar = new double[prelim_zero][prelim_one][prelim_two][prelim_three];
// 
//   int flat_count = 0;
//   
//   vector<vector<vector<vector<double> > > > recon_ar;
// 
//   for (int i = 0; i < prelim_zero; ++i) {
//     recon_ar.push_back(vector<vector<vector<double> > >());
//     for (int j = 0; j < prelim_one; ++j) {
//       recon_ar[i].push_back(vector<vector<double> >());
//       for (int k = 0; k < prelim_two; ++k) {
//         recon_ar[i][j].push_back(vector<double>());
//         for (int q = 0; q < prelim_three; ++q) {
//           recon_ar[i][j][k].push_back(0);
//         }
//       }
//     }
//   }
// 
//   for (int i = 0; i < prelim_zero; ++i) {
//     for (int j = 0; j < prelim_one; ++j) {
//       for (int k = 0; k < prelim_two; ++k) {
//         for (int q = 0; q < prelim_three; ++q) {
//           recon_ar[i][j][k][q] = flat_ar[flat_count];
//           ++flat_count;
//         }
//       }
//     }
//   }
// 
// 
//   --ax;
//   int my_dim = g_dims[ax];
//   g_dims[ax] = want_len;
// 
//   int nprelim_zero = g_dims[0];
//   int nprelim_one = g_dims[1];
//   int nprelim_two = g_dims[2];
//   int nprelim_three = g_dims[3];
// 
//   //int new_ar[nprelim_zero][nprelim_one][nprelim_two][nprelim_three];
//   vector<vector<vector<vector<double> > > > new_ar;
//   
//   for (int i = 0; i < nprelim_zero; ++i) {
//     new_ar.push_back(vector<vector<vector<double> > >());
//     for (int j = 0; j < nprelim_one; ++j) {
//       new_ar[i].push_back(vector<vector<double> >());
//       for (int k = 0; k < nprelim_two; ++k) {
//         new_ar[i][j].push_back(vector<double>());
//         for (int q = 0; q < nprelim_three; ++q) {
//           new_ar[i][j][k].push_back(0);
//         }
//       }
//     }
//   }
//   
// 
//   if (ax == 0) {
//     for (int i = 0; i < nprelim_zero; ++i) {
//       for (int j = 0; j < nprelim_one; ++j) {
//         for (int k = 0; k < nprelim_two; ++k) {
//           for (int q = 0; q < nprelim_three; ++q) {
//             new_ar[i][j][k][q] = recon_ar[0][j][k][q];
//           }
//         }
//       }
//     }
//   } else if (ax == 1) {
//     for (int i = 0; i < nprelim_zero; ++i) {
//       for (int j = 0; j < nprelim_one; ++j) {
//         for (int k = 0; k < nprelim_two; ++k) {
//           for (int q = 0; q < nprelim_three; ++q) {
//             new_ar[i][j][k][q] = recon_ar[i][0][k][q];
//           }
//         }
//       }
//     }
//   } else if (ax == 2) {
//     for (int i = 0; i < nprelim_zero; ++i) {
//       for (int j = 0; j < nprelim_one; ++j) {
//         for (int k = 0; k < nprelim_two; ++k) {
//           for (int q = 0; q < nprelim_three; ++q) {
//             new_ar[i][j][k][q] = recon_ar[i][j][0][q];
//           }
//         }
//       }
//     }
//   } else if (ax == 3) {
//     for (int i = 0; i < nprelim_zero; ++i) {
//       for (int j = 0; j < nprelim_one; ++j) {
//         for (int k = 0; k < nprelim_two; ++k) {
//           for (int q = 0; q < nprelim_three; ++q) {
//             new_ar[i][j][k][q] = recon_ar[i][j][k][0];
//           }
//         }
//       }
//     }
//   }
// 
//   int nflat_count = 0;
//   NumericVector nflat_ar(flat_count);
// 
//   for (int i = 0; i < nprelim_zero; ++i) {
//     for (int j = 0; j < nprelim_one; ++j) {
//       for (int k = 0; k < nprelim_two; ++k) {
//         for (int q = 0; q < nprelim_three; ++q) {
//           nflat_ar[nflat_count] = new_ar[i][j][k][q];
//           ++nflat_count;
//         }
//       }
//     }
//   }
//   
//   //return nflat_ar;
//   return 2;
// }
// 
// // [[Rcpp::export]]
// double receive(NumericVector lol) {
//   return lol[0];
// }



