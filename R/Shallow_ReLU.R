#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

#' Constructs the sparse approximation of h(T_i * F_i)
#' 
#' @importFrom data.table fread data.table as.data.table
#' @importFrom dplyr bind_cols
#' 
#' @export
conv_create <- function(N, K, numbase) {
  #conv = as.data.table(matrix(0L, nrow = N, ncol = K))
  
  poss = (numbase * 4) - 2
  feat <- matrix(runif(poss*K), nrow = poss, ncol = K)
  
  # Call the mapping function
  ind_map <- conv_bimap(N, numbase)
  
  conv_list <- as.data.table(feat[ind_map,])
  nb_frag = dim(conv_list)[1] / numbase
  
  # Split list then colsum them using mapping
  split_conv_list <- split(conv_list,rep(1:nb_frag, each = numbase))
  pro_conv_list <- lapply(split_conv_list, colSums)
  
  recon_conv <- bind_cols(pro_conv_list)
  recon_conv <- as.data.table(t(recon_conv))
  
  recon_conv[recon_conv < 0] <- 0
  
  return(recon_conv)
}

#' Recycles a 4D array
#' 
#' @export
four_recycle <- function(ar, ax, want_l) {

  dims <- dim(ar)
  my_dim <- dims[ax]
  dims[ax] <- want_l

  new_ar <- array(0, dim = dims)

  if (ax == 1) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
          for (q in 1:dims[4]) {
            new_ar[i, j, k, q] = ar[1, j, k, q]
          }
        }
      }
    }
  } else if (ax == 2) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
          for (q in 1:dims[4]) {
            new_ar[i, j, k, q] = ar[i, 1, k, q]
          }
        }
      }
    }
  } else if (ax == 3) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
          for (q in 1:dims[4]) {
            new_ar[i, j, k, q] = ar[i, j, 1, q]
          }
        }
      }
    }
  } else if (ax == 4) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
          for (q in 1:dims[4]) {
            new_ar[i, j, k, q] = ar[i, j, k, 1]
          }
        }
      }
    }
  }

  return (new_ar)
}





