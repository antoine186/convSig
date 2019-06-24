#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
StringMatrix simu_prep(CharacterMatrix ref, int tot_len, int omit) {
  
  int nb_line = ref.nrow();
  
  int nb_chrom = 0;
  int local_count = 0;
  int global_count = 0;
  int all_count = 0;
  int cur_line_len;
  string temp;
  string cur_char;
  //CharacterMatrix res_mat(tot_len - omti * 60, 3);
  StringMatrix res_mat(tot_len - 60, 3);
  
  for (int i = 0; i < nb_line; ++i) {
    
    temp = ref[i];
    char cur_line[temp.size() + 1];
    strncpy(cur_line, temp.c_str(), sizeof(cur_line));
    cur_line[sizeof(cur_line) - 1] = 0;
    
    cur_line_len = temp.length();
    
    if (cur_line[0] == '>') {
      ++nb_chrom;
      local_count = 0;
      continue;
    }
    
    for (int j = 0; j < cur_line_len; ++j) {
      cur_char = cur_line[j];
      
      res_mat(global_count, 0) = cur_char;
      res_mat(global_count, 1) = local_count + 1; 
      res_mat(global_count, 2) = nb_chrom; 
      
      ++global_count;
      ++local_count;
    }
    
  }
  
  return res_mat;
}
