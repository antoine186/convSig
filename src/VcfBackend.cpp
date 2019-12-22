#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
DataFrame icgc_creater(DataFrame vcf_data, NumericMatrix sample_data,
                       CharacterVector sample_names, int height, int old_height) {
  
  CharacterVector sample_new(height);
  NumericVector chrom_new(height);
  NumericVector chrom_start_new(height);
  CharacterVector ref_new(height);
  CharacterVector alt_new(height);
  
  NumericVector chrom = vcf_data[0];
  NumericVector pos = vcf_data[1];
  CharacterVector ref = vcf_data[2];
  CharacterVector alt = vcf_data[3];
  
  int df_count = 0;
  int cur_val = 0;
  int end_it = sample_data.ncol();
  
  for (int i = 0; i < end_it; i = i + 1) {

    for (int j = 0; j < old_height; j = j + 1) {

      cur_val = sample_data(j, i);
      if (cur_val == 1) {

        sample_new[df_count] = sample_names[i];
        chrom_new[df_count] = chrom[j];
        chrom_start_new[df_count] = pos[j];
        ref_new[df_count] = ref[j];
        alt_new[df_count] = alt[j];

        ++df_count;

      }

    }

  }
  
  DataFrame df = DataFrame::create( Named("icgc_sample_id") = sample_new,
                                    Named("chromosome") = chrom_new,
                                    Named("chromosome_start") = chrom_start_new,
                                    Named("chromosome_end") = chrom_start_new,
                                    Named("mutated_from_allele") = ref_new,
                                    Named("mutated_to_allele") = alt_new);
  
  return df;
  
}

