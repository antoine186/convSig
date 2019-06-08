#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

#' Constructs the sparse approximation of h(T_i * F_i)
#' 
#' @importFrom data.table fread data.table as.data.table
#' @importFrom dplyr bind_cols
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
  
  invisible(recon_conv)
}

# Recycles a 4D array at the index specified by the user
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

  invisible(new_ar)
}

# Stores the indices of the fragment types according to which base is present
# at a particular base position
fragbase_indexer <- function(numbase, N) {
  
  type <- list()
  base_div = 1
  
  if (numbase == 3) {
    mid = 2
  } else if (numbase == 5) {
    mid = 3
  }

  for(i in 1:numbase) {
    
    if (i == mid) {
      temp_struct <- list(c(), c())
      base_mod = 2
    } else {
      temp_struct <- list(c(), c(), c(), c())
      base_mod = 4
    }
    
    for (j in 1:N) {
      inter_j = j - 1
      
      base_type = (floor(inter_j / base_div) %% base_mod) + 1
      temp_struct[[base_type]] <- c(temp_struct[[base_type]], j) 
    }
    
    type[[i]] <- temp_struct
    
    if (i == mid) {
      base_div <- base_div * 2
    } else {
      base_div <- base_div * 4
    }
    
  }
  
  invisible(type)
}

# Returns a matrix with all possible unmutated fragment types one-hot encoded
tencode <- function(numbase, N) {
  
  T = matrix(0, nrow = N, ncol = (4 * numbase) - 2)
  
  if (numbase == 3) {
    mid = 2
  } else if (numbase == 5) {
    mid = 3
  }
  
  for (i in 1:N) {
    #temp = i
    inter_temp = i - 1
    for (j in 1:numbase) {
      inter_j = j - 1
      
      if (j < mid) {
        T[i,(4 * inter_j + inter_temp %% 4) + 1]=1
        
        #T[i,4*j+temp%%4]=1;
        inter_temp= floor(inter_temp / 4)
      } else if (j == mid) {
        T[i,(4 * inter_j + inter_temp %% 2) + 1]=1
        inter_temp = floor(inter_temp / 2)
      } else if (j > mid) {
        T[i,(4 * inter_j - 2 + inter_temp %% 4) + 1]=1
        inter_temp = inter_temp / 4
      }
      
    }
  }
  
  invisible(T)
}

# Imports and processes the mutation count data
mutation_inputprocess <- function(X, numbase) {

  S <- dim(X)[1]
  N = 4^(numbase - 1) * 2

  X <- c(t(as.matrix(X)))
  X_ar <- array(0, dim = c(S, N, 3))

  count = 1
  for (i in 1:S) {
    for (j in 1:N) {
      for (k in 1:3) {
        X_ar[i,j,k] = X[count]
        count = count + 1
      }
    }
  }

  invisible(X_ar)
}

# Imports and processes the mutation count data
background_inputprocess <- function(X, numbase) {
  
  N = 4^(numbase - 1) * 2
  
  X <- c(t(as.matrix(X)))
  X_ar <- array(0, dim = c(N, 3))
  
  count = 1
  for (i in 1:N) {
    for (j in 1:3) {
      X_ar[i,j] = X[count]
      count = count + 1
    }
  }
  X_ar <- X_ar[,1]
  
  invisible(X_ar)
}

# Splits the counts found in the mutation count file into test/train counts
test_splitX <- function(X_input, S, N) {
  
  X_ar <- array(0, dim = c(S, N, 3))
  for (i in 1:S) {
    for (j in 1:N) {
      for (k in 1:3) {
        X_ar[i,j,k] = Xsplitter(X_input[i,j,k])
      }
    }
  }
  
  invisible(X_ar)
}

# Splits the counts found in the background count file into test/train counts
test_splitbg <- function(bg_input, numbase) {
  
  if (numbase == 3) {
    type_length = 32
  } else if (numbase == 5) {
    type_length = 512
  }
  
  X_ar <- array(0, dim = type_length)
  for (i in 1:type_length) {
    X_ar[i] = Xsplitter(bg_input[i])
  }
  
  invisible(X_ar)
}

Xsplitter <- function(count_nb) {
  
  test_X <- rbinom(1, count_nb, 0.5)
  
  invisible(test_X)
}

# Sums a 4D matrix at the 4th axis
four_colsum <- function(X, ax) {
  
  dims <- dim(X)
  
  i_final <- dim(X)[1]
  j_final <- dim(X)[2]
  k_final <- dim(X)[3]
  q_final <- dim(X)[4]
  
  dims[ax] = 1
  
  X_ar <- array(0, dim = dims)
  
  if (ax == 4) {
    for (i in 1:i_final) {
      for (j in 1:j_final) {
        for (k in 1:k_final) {
          acc <- array(0, dim = q_final)
          for (q in 1:q_final) {
            acc[q] = X[i,j,k,q]
          }
          
          X_ar[i,j,k,1] = sum(acc)
        }
      }
    }
  }
  
  invisible(X_ar)
}

# Sums a 3D matrix at the 2nd axis
three_colsum <- function(X, ax) {
  
  dims <- dim(X)
  
  i_final <- dim(X)[1]
  j_final <- dim(X)[2]
  k_final <- dim(X)[3]
  
  dims[ax] = 1
  
  X_ar <- array(0, dim = dims)
  
  if (ax == 2) {
    for (i in 1:i_final) {
      for (k in 1:k_final) {
        acc <- array(0, dim = j_final)
        for (j in 1:j_final) {
          acc[j] = X[i,j,k] 
        }
        
        X_ar[i,1,k] = sum(acc)
      }
    }
  }
  
  invisible(X_ar)
}

# Function that creates the Z hidden variable for the EM algorithm
hidden_create <- function(S, N, K) {
  
  Z <- array(runif(S*N*3*K), dim = c(S, N, 3, K))
  Zk_sum <- four_colsum(Z, 4)
  
  res <- sweep(Z,MARGIN=c(1,2,3),Zk_sum,`/`)
  
  invisible(res)
}

#' Solves the |T * Feature - theta|^2 = 0 linear regression problem
#' 
#' @importFrom MASS ginv
linreg_solver <- function(T, theta) {
  
  T_inv <- ginv((t(T) %*% T))
  Multi_F = T_inv %*% t(T)
  feat = Multi_F %*% theta
  
  invisible(feat)
}

#' Performs the ReLU transform on the mutational data
#' 
#' @export
relu_transform <- function(mut_obj, five = FALSE, K = 5) {
  
  if (K == 0) {
    stop("Your specified number of mutational processes cannot be zero")
  } else if (K > 80) {
    stop("Your specified number of mutational processes cannot exceed 80")
  }
  
  if (!("Shallowres" %in% class(mut_obj))) {
    stop("Your specified mut_obj is not of the right type")
  }
  
  if (five == FALSE) {
    numbase = 3
  } else if (five == TRUE) {
    numbase = 5
  }
  
  mut_path <- mut_obj@mut_mat
  bg_path <- mut_obj@wt
  
  X <- mutation_inputprocess(mut_path, numbase)
  bg <- background_inputprocess(bg_path, numbase)
  
  S <- dim(X)[1]
  N <- dim(X)[2]
  
  X_test <- test_splitX(X, S, N)
  bg_test <- test_splitbg(bg, numbase)
  
  Z <- hidden_create(S, N, K)
  
  P <- array(runif(S*K), dim = c(S, K))
  
  conv <- conv_create(N, K, numbase)
  
  theta_array = array(runif(N*K), dim = c(N, K))
  theta_array <- sweep(theta_array,MARGIN=c(2),colSums(theta_array),`/`)
  
  beta_array = X + 10^(-4)
  
  mat <- array(runif(2*3*K), dim = c(2, 3, K))
  summed_mat <- three_colsum(mat, 2)
  mat <- sweep(mat,MARGIN=c(1,3),summed_mat, `/`)
  
  type <- fragbase_indexer(numbase, N)
  T <- tencode(numbase, N)
  
  
  
  invisible(mat)
}

# initialize the initial LOSS value
init_LOSS <- function(bg, theta, P) {
  
  bg <- array(bg, dim = c(length(bg), 1))
  LOSS <- sweep(theta,MARGIN=c(1),bg, `*`)
  
  invisible(LOSS) 
}

# Loop function, which increments lambda (burden of the loss function) iteratively
regularizer <- function() {
  
  for (r in 1:6) {
    
    reg = (10^r) - 1
    
    
    
  }
  
}


