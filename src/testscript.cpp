#include <Rcpp.h>
#include <string>
#include <cstring>

using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
RcppExport SEXP timesTwo(SEXP x) {
  double res = as<double>(x) * 2;
  throw "lol";
  return wrap(res);
}

// [[Rcpp::export]]
RcppExport SEXP timesTwoList(SEXP x) {
  NumericVector res = as<NumericVector>(x) * 2;
  return wrap(res);
}

// [[Rcpp::export]]
NumericVector WeirdVector(NumericVector x) {
  std::vector<int> temp_array;
  temp_array.push_back(3); 
  temp_array.push_back(4);
  temp_array.push_back(5);
  
  x[temp_array[0]] = ++x[temp_array[0]];
  x[temp_array[1]] = ++x[temp_array[1]];
  return x;
}

// [[Rcpp::export]]
bool string_comptest(CharacterVector x) {
  bool yo;
  if (x[0] == "hello") {
    yo = true;
  } else {
    yo = false;
  }
  return yo;
}

// [[Rcpp::export]]
int string_comptest2(CharacterVector x, CharacterVector y, int a) {
  int index = 0;
  for(CharacterVector::iterator it = x.begin(); 
      it != x.end(); ++it) {
    //++index;
    if (*it == y[a]) {
    //if (strncmp(it.c_str(), "lol") == 1) {
    //if (x.find("lol") != std::string::npos) {
      //++index;
      index = it - x.begin();
      break;
    }
  }
  return index;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
