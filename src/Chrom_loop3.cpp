#include <Rcpp.h>
#include <string>
#include <regex>
#include <cstring>
#include <array>
#include <unordered_map>
using namespace Rcpp;
using namespace std;

// This function processes the fasta file along with the mutation file

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
S4 shallow_loop3(S4 mat) {
  
  NumericMatrix mut_mat = mat.slot("mut_mat");
  mut_mat(0,1) = 5;
  mut_mat(0,2) = 4;
  mut_mat(0,3) = 3;
  mat.slot("mut_mat") = mut_mat;
  //mat.slot("name")  = "Sewall Wright";
  return mat;
}

