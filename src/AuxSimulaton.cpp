#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
String simu_prep(CharacterMatrix ref) {
  
  int nb_line = ref.nrow();
  
  int nb_chrom = 0;
  int cur_line_len;
  string temp;
  //string cur_line;
  
  for (int i = 0; i < nb_line; ++i) {
    
    temp = ref[i];
    char cur_line[temp.size() + 1];
    strncpy(cur_line, temp.c_str(), sizeof(cur_line));
    cur_line[sizeof(cur_line) - 1] = 0;
    
    cur_line_len = temp.length();
    
    if (cur_line[0] == '>') {
      ++nb_chrom;
      continue;
    } 
    
  }
  
  return nb_chrom;
}

// else {
//   for (int j = 0; j <cur_line_len; ++j) {
//     cur_line[j]
//   }
// }
