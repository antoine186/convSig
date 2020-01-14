#include <Rcpp.h>
#include <string>
#include <regex>
#include <cstring>
#include <array>
#include <unordered_map>
#include <typeinfo>
using namespace Rcpp;
using namespace std;

// This function processes the fasta file along with the mutation file

// Base to number function
int base2num(std::string base) {
  if (base == "A") {
    return(0);
  } else if (base == "G") {
    return(1);
  } else if (base == "C") {
    return(2);
  } else if (base == "T") {
    return(3);
  } else {
    return(4);
  }
}

// [[Rcpp::export]]
std::vector<int> feat2table3(std::string b1, std::string b2, std::string b3) {
  int n1 = base2num(b1);
  int n2 = base2num(b2);
  int n3 = base2num(b3);
  int n1t = n1;
  
  if (n2 < 2) {
    n1 = 3 - n3;
    n3 = 3 - n1t;
  } else {
    n2 = 3 - n2;
  }
  
  int value = (n1*8*3)+(n2*4*3)+(n3*3);
  vector<int> values;  
  values.push_back(value);
  values.push_back(value + 1);
  values.push_back(value + 2);
  
  return(values);
}

// [[Rcpp::export]]
std::vector<int> feat2table5(std::string b1, std::string b2, std::string b3,
                             std::string b4, std::string b5) {
  
  int n1 = base2num(b1);
  int n2 = base2num(b2);
  int n3 = base2num(b3);
  int n4 = base2num(b4);
  int n5 = base2num(b5);
  
  int n1t = n1;
  int n2t = n2;
  
  if (n3 < 2) {
    n1 = 3 - n5;
    n2 = 3 - n4;
    n4 = 3 - n2t;
    n5 = 3 - n1t;
  } else {
    n3 = 3 - n3;
  }
  
  int value = (n1*128*3)+(n2*32*3)+(n3*16*3)+(n4*4*3)+(n5*3);
  vector<int> values;  
  values.push_back(value);
  values.push_back(value + 1);
  values.push_back(value + 2);
  
  return(values);
}

// [[Rcpp::export]]
DataFrame reverse_transform(DataFrame alleles) {
  NumericVector from_allele = alleles[0];
  NumericVector to_allele = alleles[1];
  
  int nr = alleles.nrow();
  
  for (int i = 0; i < nr; ++i)  
  { 
    if (from_allele[i] < 2) {
      to_allele[i] = 3 - to_allele[i];
    }
    
    if (to_allele[i] > 2) {
      to_allele[i] = 2;
    }
  }
  
  DataFrame mod_alleles = DataFrame::create( Named("from_allele") = from_allele,  
                                             Named("to_allele") = to_allele);
  
  return mod_alleles;
}

// [[Rcpp::export]]
S4 shallow_loop3(S4 mat, DataFrame fasta, DataFrame mut_file, CharacterVector uniq_samples) {
  
  NumericMatrix mut_mat = mat.slot("mut_mat");
  NumericVector wt_ar = mat.slot("wt");
  
  CharacterVector sampleid_ar = mut_file[0];
  NumericVector chromid_ar = mut_file[1];
  NumericVector startpos_ar = mut_file[2];
  CharacterVector ref_ar = mut_file[3];
  NumericVector mut_ar = mut_file[4];
  
  int loop1 = fasta.nrow();
  int mut_file_length = mut_file.nrow();
  
  CharacterVector f = as<CharacterVector>(fasta[0]);
  
  int chrom = 0;
  int ref_pos = -1;
  int mut_pos = 0;
  std::string base1;
  std::string base2;
  std::string base3;
  
  //int fasta_stat = 0;
  //int spe_stat = 0;
  //int line_stat = 0;
  //int feature_stat = 0;
  //int free_stat = 0;
  
  std::regex chrompattern ("^>");
  // For loop 1
  for (int i = 0; i < loop1; i = i + 1) {
    std::string cur_line = as<std::string>(f[i]);
    //++fasta_stat;
    // Special case for chromosome header line
    if (std::regex_search(cur_line, chrompattern)) {
      chrom++;
      
      if (chrom > 24) {
        break;
      }
      
      ref_pos = -1;
      base1 = "N";
      base2 = "N";
      continue;
    }
    
    // Normal case for line with nucleotides
    int n = cur_line.length(); 
    char cur_array[n + 1]; 
    strcpy(cur_array, cur_line.c_str());
    
    // For loop 2 based on cur_array
    int loop2 = strlen(cur_array); 
    for (int j = 0; j < loop2; j = j + 1) {
      //++line_stat;
      base3 = cur_array[j];
      std::vector<int> temp_array;
      
      if (base1 != "N" && (base2 != "N" && base3 != "N")) {
        
        temp_array = feat2table3(base1, base2, base3);
        
        ++wt_ar[temp_array[0]];
        ++wt_ar[temp_array[1]];
        ++wt_ar[temp_array[2]];
        
        //++feature_stat;
      }
      
      ++ref_pos;
      
      if (ref_pos > mut_pos) {
        std::string err_ref = std::to_string(ref_pos);
        std::string err_mut = std::to_string(mut_pos);
        
        stop(err_ref + " " + err_mut);
      }
      
      // Here is the second while loop - startpos_ar[mut_pos] is generally
      // larger than ref_pos
      if (mut_pos < mut_file_length && !(chromid_ar[mut_pos] > chrom ||
      startpos_ar[mut_pos] > ref_pos)) {
        //++free_stat;
        while (1) {
          //++free_stat;
          if(ref_ar[mut_pos] != base2) {
            stop("It seems that your supplied assembly input does not "
                   "correspond to the one in your mutation input file");
          }
          
          int index;
          
          for(CharacterVector::iterator it = uniq_samples.begin(); 
              it != uniq_samples.end(); ++it) {
            
            if (*it == sampleid_ar[mut_pos]) {
              //++free_stat;
              index = it - uniq_samples.begin();
              break;
            }
            
          }
          
          ++mut_mat(index, temp_array[mut_ar[mut_pos]]);
          
          mut_pos = ++mut_pos;
          
          if (mut_pos >= mut_file_length || (chromid_ar[mut_pos] > chrom ||
              startpos_ar[mut_pos] > ref_pos)) {
            break;
          }
        }
      } 
      
      base1 = base2;
      base2 = base3;
      
      std::vector<int>().swap(temp_array);
      
    }
  }
  
  mat.slot("mut_mat") = mut_mat;
  mat.slot("wt") = wt_ar;
  //mat.slot("fasta_status") = fasta_stat;
  //mat.slot("special_status") = spe_stat;
  //mat.slot("perline_status") = line_stat;
  //mat.slot("feature_status") = feature_stat;
  //mat.slot("free_status") = free_stat;
  
  return mat;
  
}

// [[Rcpp::export]]
S4 shallow_loop5(S4 mat, DataFrame fasta, DataFrame mut_file, CharacterVector uniq_samples) {
  
  NumericMatrix mut_mat = mat.slot("mut_mat");
  NumericVector wt_ar = mat.slot("wt");
  
  CharacterVector sampleid_ar = mut_file[0];
  NumericVector chromid_ar = mut_file[1];
  NumericVector startpos_ar = mut_file[2];
  CharacterVector ref_ar = mut_file[3];
  NumericVector mut_ar = mut_file[4];
  
  int loop1 = fasta.nrow();
  int mut_file_length = mut_file.nrow();
  
  CharacterVector f = as<CharacterVector>(fasta[0]);
  
  int chrom = 0;
  int ref_pos = -2;
  int mut_pos = 0;
  std::string base1;
  std::string base2;
  std::string base3;
  std::string base4;
  std::string base5;
  
  //int fasta_stat = 0;
  //int spe_stat = 0;
  //int line_stat = 0;
  //int feature_stat = 0;
  //int free_stat = 0;
  
  std::regex chrompattern ("^>");
  // For loop 1
  for (int i = 0; i < loop1; i = i + 1) {
    std::string cur_line = as<std::string>(f[i]);
    //++fasta_stat;
    // Special case for chromosome header line
    if (std::regex_search(cur_line, chrompattern)) {
      chrom++;
      
      if (chrom > 24) {
        break;
      }
      
      ref_pos = -2;
      base1 = "N";
      base2 = "N";
      base3 = "N";
      base4 = "N";
      continue;
    }
    
    // Normal case for line with nucleotides
    int n = cur_line.length(); 
    char cur_array[n + 1]; 
    strcpy(cur_array, cur_line.c_str());
    
    // For loop 2 based on cur_array
    int loop2 = strlen(cur_array); 
    for (int j = 0; j < loop2; j = j + 1) {
      //++line_stat;
      base5 = cur_array[j];
      std::vector<int> temp_array;
      
      if (base1 != "N" && (base2 != "N" && (base3 != "N" && (base4 != "N" && base5 != "N")))) {
        
        temp_array = feat2table5(base1, base2, base3, base4, base5);
        
        ++wt_ar[temp_array[0]];
        ++wt_ar[temp_array[1]];
        ++wt_ar[temp_array[2]];
        
        //++feature_stat;
      }
      
      ++ref_pos;
      
      // Here is the second while loop - startpos_ar[mut_pos] is generally
      // larger than ref_pos
      if (mut_pos < mut_file_length && !(chromid_ar[mut_pos] > chrom ||
      startpos_ar[mut_pos] > ref_pos)) {
        //++free_stat;
        while (1) {
          //++free_stat;
          if(ref_ar[mut_pos] != base3) {
            stop("It seems that your supplied assembly input does not "
                   "correspond to the one in your mutation input file");
          }
          
          int index;
          
          for(CharacterVector::iterator it = uniq_samples.begin(); 
              it != uniq_samples.end(); ++it) {
            
            if (*it == sampleid_ar[mut_pos]) {
              //++free_stat;
              index = it - uniq_samples.begin();
              break;
            }
            
          }
          
          ++mut_mat(index, temp_array[mut_ar[mut_pos]]);
          
          mut_pos = ++mut_pos;
          
          if (mut_pos >= mut_file_length || (chromid_ar[mut_pos] > chrom ||
              startpos_ar[mut_pos] > ref_pos)) {
            break;
          }
        }
      } 
      
      base1 = base2;
      base2 = base3;
      base3 = base4;
      base4 = base5;
      
      std::vector<int>().swap(temp_array);
      
    }
  }
  
  mat.slot("mut_mat") = mut_mat;
  mat.slot("wt") = wt_ar;
  //mat.slot("fasta_status") = fasta_stat;
  //mat.slot("special_status") = spe_stat;
  //mat.slot("perline_status") = line_stat;
  //mat.slot("feature_status") = feature_stat;
  //mat.slot("free_status") = free_stat;
  
  return mat;
  
}