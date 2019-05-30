#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

#' Constructs the sparse approximation of h(T_i * F_i)
#' 
#' @importFrom data.table fread data.table as.data.table
conv_create <- function(N, K, numbase) {
  conv = as.data.table(matrix(0L, nrow = N, ncol = K))
  
  feat <- matrix(runif(18*5), nrow = 18, ncol = 5)
  
  # Call the mapping function
  ind_map <- conv_bimap
}