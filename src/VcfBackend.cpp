#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
DataFrame icgc_creater(DataFrame vcf_data, CharacterVector cnames, int height, 
                    int start_it, int end_it) {
  
  NumericVector num_vec(height);
  CharacterVector char_vec(height);

  // When giving names to columns
  DataFrame df = DataFrame::create( Named("icgc_sample_id") = char_vec , 
                                    Named("chromosome") = num_vec, 
                                    Named("chromosome_start") = num_vec, 
                                    Named("chromosome_end") = num_vec, 
                                    Named("mutated_from_allele") = char_vec, 
                                    Named("mutated_to_allele") = char_vec);
  
  int df_count = 0;
  int cur_val = 0;
  
  for (int i = start_it; i < end_it + 1; i = i + 1) {

    for (int j = 0; j < height; j = j + 1) {
      
      cur_val = vcf_data[j, i];
      if (cur_val == 1) {

        df[df_count, 0] = cnames[i];
        df[df_count, 1] = vcf_data[j, 0];
        df[df_count, 2] = vcf_data[j, 1];
        df[df_count, 4] = vcf_data[j, 3];
        df[df_count, 5] = vcf_data[j, 4];
        
        ++df_count;

      }

    }

  }
  
  return df;
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
