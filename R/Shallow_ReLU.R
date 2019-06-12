#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

setClass (
  # Class name
  "feat_obj",
  
  # Defining slot type
  representation (
    feat = "matrix",
    Multi_F = "matrix"
  )
)

#' Constructs the sparse approximation of h(T_i * F_i)
#' 
#' @importFrom data.table fread data.table as.data.table
#' @importFrom dplyr bind_cols
conv_create <- function(N, K, numbase, feat) {
  #conv = as.data.table(matrix(0L, nrow = N, ncol = K))
  
  # poss = (numbase * 4) - 2
  # feat <- matrix(runif(poss*K), nrow = poss, ncol = K)
  
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
  
  invisible(as.matrix(recon_conv))
}

# Recycles a 3D array at the index specified by the user
three_recycle <- function(ar, ax, want_l) {
  
  dims <- dim(ar)
  my_dim <- dims[ax]
  dims[ax] <- want_l
  
  new_ar <- array(0, dim = dims)
  
  if (ax == 1) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
          new_ar[i, j, k] = ar[1, j, k]
        }
      }
    }
  } else if (ax == 2) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
            new_ar[i, j, k] = ar[i, 1, k]
        }
      }
    }
  } else if (ax == 3) {
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
            new_ar[i, j, k] = ar[i, j, 1]
        }
      }
    }
  }
  
  invisible(new_ar)
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

# Sums a 4D matrix at the 2nd, 3rd and 4th axis
four_colsum <- function(X, ax) {
  
  dims <- dim(X)
  
  i_final <- dim(X)[1]
  j_final <- dim(X)[2]
  k_final <- dim(X)[3]
  q_final <- dim(X)[4]
  
  dims[ax] = 1
  
  X_ar <- array(0, dim = dims)
  
  if (ax == 2) {
    for (i in 1:i_final) {
      for (k in 1:k_final) {
        for (q in 1:q_final) {
          acc <- array(0, dim = j_final)
          for (j in 1:j_final) {
            acc[j] = X[i,j,k,q]
          }
          
          X_ar[i,1,k,q] = sum(acc)
        }
      }
    }
  } else if (ax == 3) {
    for (i in 1:i_final) {
      for (j in 1:j_final) {
        for (q in 1:q_final) {
          acc <- array(0, dim = k_final)
          for (k in 1:k_final) {
            acc[k] = X[i,j,k,q]
          }
          
          X_ar[i,j,1,q] = sum(acc)
        }
      }
    }
  } else if (ax == 4) {
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

# Sums a 5D matrix at the 3rd axis
five_colsum <- function(X, ax) {
  
  dims <- dim(X)
  
  i_final <- dim(X)[1]
  j_final <- dim(X)[2]
  k_final <- dim(X)[3]
  q_final <- dim(X)[4]
  y_final <- dim(X)[5]
  
  dims[ax] = 1
  
  X_ar <- array(0, dim = dims)
  
  if (ax == 3) {
    for (i in 1:i_final) {
      for (j in 1:j_final) {
        for (y in 1:y_final) {
          for (q in 1:q_final) {
            acc <- array(0, dim = k_final)
            for (k in 1:k_final) {
              acc[k] = X[i,j,k,q,y]
            }
            
            X_ar[i,j,1,q,y] = sum(acc)
          }
        }
      }
    }
  }
  
  invisible(X_ar)
}

# Sums a 3D matrix at the 2nd or 3rd axis
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
  } else if (ax == 3) {
    for (i in 1:i_final) {
      for (j in 1:j_final) {
        acc <- array(0, dim = k_final)
        for (k in 1:k_final) {
          acc[k] = X[i,j,k] 
        }
        
        X_ar[i,j,1] = sum(acc)
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
  
  feat_o <- new("feat_obj", feat = feat, Multi_F = Multi_F)
  
  invisible(feat_o)
}

# A function which helps indexing a 3D matrix using complex list indexing
list_indexer <- function(ar, type, mid, N, K, X = TRUE) {
  
  dims <- dim(ar)
  list1 <- unlist(type[[mid]][1])
  list2 <- unlist(type[[mid]][2])
  
  mut_nb <- N/2
  
  if (X == TRUE) {
    new_ar <- array(0, dim = c(dims[1], 2, mut_nb, 3))
    
    for (i in 1:dims[1]) {
      for (j in 1:2) {
        list_mut <- unlist(type[[mid]][j])
        for (k in 1:mut_nb) {
          fancy_ind = list_mut[k]
          for (q in 1:3) {
            new_ar[i,j,k,q] = ar[i, fancy_ind, q]
          }
        }
      }
    }
  } else {
    new_ar <- array(0, dim = c(dims[1], 2, mut_nb, 3, K))
    
    for (i in 1:dims[1]) {
      for (j in 1:2) {
        list_mut <- unlist(type[[mid]][j])
        for (k in 1:mut_nb) {
          fancy_ind = list_mut[k]
          for (q in 1:3) {
            for (y in 1:K) {
              new_ar[i,j,k,q,y] = ar[i, fancy_ind, q, y]
            }
          }
        }
      }
    }
  }
  
  invisible(new_ar)
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
    mid = 2
  } else if (five == TRUE) {
    numbase = 5
    mid = 3
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
  
  poss = (numbase * 4) - 2
  feat <- matrix(runif(poss*K), nrow = poss, ncol = K)
  
  conv <- conv_create(N, K, numbase, feat)
  
  theta_array = array(runif(N*K), dim = c(N, K))
  theta_array <- sweep(theta_array,MARGIN=c(2),colSums(theta_array),`/`)
  
  beta_array = X + 10^(-4)
  
  mat <- array(runif(2*3*K), dim = c(2, 3, K))
  summed_mat <- three_colsum(mat, 2)
  mat <- sweep(mat,MARGIN=c(1,3),summed_mat, `/`)
  
  type <- fragbase_indexer(numbase, N)
  T <- tencode(numbase, N)
  
  feat_o <- linreg_solver(T, theta_array)
  
  reg_res <- regularizer(X, bg, conv, theta_array, P, mat, N, S, K, type,
                         mid, beta_array, Z,
                         feat_o@Multi_F, feat_o@feat, numbase, bg_test, X_test, T)
  
  invisible(reg_res)
}

# initialize the initial LOSS value
init_LOSS <- function(bg, theta, P) {
  
  bg <- array(bg, dim = c(length(bg), 1))
  inter_LOSS <- sweep(theta,MARGIN=c(1),bg, `*`)
  inter_LOSS <- sweep(P,MARGIN=c(2),colSums(inter_LOSS),`*`)
  LOSS <- sum(inter_LOSS)
  
  invisible(LOSS) 
}

# Loop function, which increments lambda (burden of the loss function) iteratively
regularizer <- function(X, bg, conv, theta, P, mat, N, S, K,
                        type, mid, beta_ar, Z, Multi_F, feat, numbase, bg_test, X_test, T) {
  
  for (r in 1:6) {
    
    reg = (10^(r - 1)) - 1
    
    cat("Using lambda value: ")
    cat(reg)
    cat("\n")
    
    LOSS <- init_LOSS(bg, theta, P)
    for (i in 1:2) {
      
      #
      inter_theta <- theta[unlist(type[[mid]][i]),]
      theta_dim <- dim(inter_theta)
      inter_theta <- array(inter_theta, dim = c(theta_dim[1],1,K))
      
      inter_P <- array(P, dim = c(S,1,1,K))
      inter_P <- four_recycle(inter_P, 2, theta_dim[1])
      
      inter_mat <- array(mat[i,,], dim = c(1,3,K))
      inter_mat <- three_recycle(inter_mat, 1, theta_dim[1])
      
      temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_theta, `*`)
      temp <- four_recycle(temp, 3, 3)
      temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
      temp <- four_colsum(temp, 4)
      #
      
      inter_bg <- array(bg[unlist(type[[mid]][i])], dim = c(1,theta_dim[1],1))
      inter_bg <- three_recycle(inter_bg, 3, 3)
      inter_bg <- three_recycle(inter_bg, 1, 10)
      
      temp <- log(sweep(temp,MARGIN=c(1,2,3),inter_bg, `*`))
      
      LOSS = LOSS - sum(sweep(temp,MARGIN=c(1,2,3),X[,unlist(type[[mid]][i]),], `*`))
      
    }
    LOSS = 0.5 * reg * sum((theta - conv)^2)
    
    old_LOSS = LOSS + 1
    new_LOSS = LOSS
    
    while(abs(new_LOSS - old_LOSS) > (1.0/(2*r+1))) {
      cat("Optimising on the loss function \n")
      for (i in 1:K) {
        
        C = reg * conv[,i] - bg * sum(P[,i])
        temp2 <- array(0, dim = c(S, N, 3))
        
        for (j in 1:2) {
          
          inter_theta <- theta[unlist(type[[mid]][j]),]
          theta_dim <- dim(inter_theta)
          inter_theta <- array(inter_theta, dim = c(theta_dim[1],1,K))
          
          inter_P <- array(P, dim = c(S,1,1,K))
          inter_P <- four_recycle(inter_P, 2, theta_dim[1])
          
          inter_mat <- array(mat[j,,], dim = c(1,3,K))
          inter_mat <- three_recycle(inter_mat, 1, theta_dim[1])
          
          temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_theta, `*`)
          temp <- four_recycle(temp, 3, 3)
          temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
          temp <- four_colsum(temp, 4)
          
          temp2[,unlist(type[[mid]][j]),] <- temp
          
        }
        alpha_array <- temp2
        
        #r
        for (j in 1:2) {
          
          inter_theta <- theta[unlist(type[[mid]][j]),i]
          theta_len <- length(inter_theta)
          inter_theta <- array(inter_theta, dim = c(theta_len,1))
          
          inter_P <- array(P[,i], dim = c(S,1,1))
          inter_P <- three_recycle(inter_P, 2, theta_len)
          
          inter_mat <- array(mat[j,,i], dim = c(3))
          
          temp <- sweep(inter_P,MARGIN=c(2,3),inter_theta, `*`)
          temp <- three_recycle(temp, 3, 3)
          temp <- sweep(temp,MARGIN=c(3),inter_mat, `*`)
          
          inter_P <- array(P[,i], dim = c(S,1,1))
          inter_P <- three_recycle(inter_P, 3, 3)
          
          temp2 <- sweep(inter_P,MARGIN=c(3),inter_mat, `*`)
          
          alpha_array[,unlist(type[[mid]][j]),] = 
            alpha_array[,unlist(type[[mid]][j]),] - temp
          alpha_array[,unlist(type[[mid]][j]),] = 
            alpha_array[,unlist(type[[mid]][j]),] / three_recycle(temp2, 2, theta_len)
          
        }
        
        # sum(0,2)
        inter_upper = three_colsum((beta_ar / alpha_array), 3)
        upper = colSums(inter_upper, 1)
        upper = upper + C
        
        v_max = max(upper)
        x_max_max = array(0, dim = c(N))
        x_max_min = array(0, dim = c(N))
        
        inter_upper = three_colsum((beta_ar / (alpha_array + 1)), 3)
        upper2 = colSums(inter_upper, 1) - reg + C
        
        v_min = max(upper2)
        
        x_min_max = array(1, dim = c(N))
        x_min_min = array(0, dim = c(N))
        x_min_min[which(upper2 == max(upper2))] = 1
        
        while ((v_max - v_min) > 0.0001) {
          
          v_new = (1 / 2) * (v_max + v_min)
          x_new_max = (1 / 2) * (x_min_max + x_max_max)
          
          if (!is.finite(v_max)) {
            
            v_new = 2 * abs(v_min)
            x_new_max = x_min_max
          
          }
          
          x_new_max[v_new > upper] = 0 
          x_new_min = x_max_min
          x_new_min[v_new > upper] = 0 
          
          while ((sum(x_new_min) - 1) * (sum(x_new_max) - 1) < 0) {
            
            x_new = (1 / 2) * (x_new_min + x_new_max)
            
            inter <- three_colsum((beta_ar / sweep(alpha_array,
                                                   MARGIN=c(2),x_new, `+`)), 3)
            b_array = (array(colSums(inter, 1), dim = N) - reg * x_new + C) > v_new
            x_new_min[b_array] = x_new[b_array]
            
            x_new_max[!b_array] = x_new[!b_array]
            
          }
          
          if (sum(x_new_min) - 1 >= 0 & sum(x_new_max) - 1 >= 0) {
            
            v_min = v_new
            x_min_min = x_new_min
            x_min_max = x_new_max
            
          } else {
            
            v_max = v_new
            x_max_min = x_new_min
            x_max_max = x_new_max
            
          }
          
        }
        
        v_new = (1 / 2) * (v_max + v_min)
        x_new_max = (1 / 2) * (x_min_max + x_max_max)
        x_new_max[v_new > upper] = 0
        x_new_min = x_max_min
        x_new_min[v_new > upper] = 0
        
        while (sum((x_new_min - x_new_max)^2) > 0.001) {
          
          x_new = 1 / 2 * (x_new_min + x_new_max)
          
          inter <- three_colsum((beta_ar / sweep(alpha_array,
                                                 MARGIN=c(2),x_new, `+`)), 3)
          b_array = (array(colSums(inter, 1), dim = N) - reg * x_new + C) > v_new
          
          x_new_min[b_array] = x_new[b_array]
          x_new_max[!b_array] = x_new[!b_array]
          
        }
        
        theta[,i] = (1 / 2) * (x_new_max + x_new_min)
        
      }
      
      theta = theta / colSums(theta)
      
      for (i in 1:2) {
        
        inter_theta <- theta[unlist(type[[mid]][i]),]
        theta_dim <- dim(inter_theta)
        inter_theta <- array(inter_theta, dim = c(theta_dim[1],1,K))
        
        inter_P <- array(P, dim = c(S,1,1,K))
        inter_P <- four_recycle(inter_P, 2, theta_dim[1])
        
        inter_mat <- array(mat[i,,], dim = c(1,3,K))
        inter_mat <- three_recycle(inter_mat, 1, theta_dim[1])
        
        temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_theta, `*`)
        temp <- four_recycle(temp, 3, 3)
        temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
        
        Z[,unlist(type[[mid]][i]),,] = temp

      }
      
      Z = Z + (1e-4 * sum(Z) / (K*N*3*S))
      Z = sweep(Z,MARGIN=c(1,2,3),four_colsum(Z, 4), `/`)
      
      inter_X = list_indexer(X, type, mid, N, K)
      inter_Z = list_indexer(Z, type, mid, N, K, X = FALSE)
        
      inter_mat = five_colsum(sweep(inter_Z,MARGIN=c(1,2,3,4),inter_X, `*`), 3)
      mat = array((colSums(inter_mat) + 0.01), dim = c(2,3,K))
      summed_mat <- three_colsum(mat, 2)
      mat <- sweep(mat,MARGIN=c(1,3),summed_mat, `/`)
      
      inter_P = sweep(Z,MARGIN=c(1,2,3),X, `*`)
      inter_P = four_colsum(inter_P, 2)
      inter_P = array(four_colsum(inter_P, 3), dim = c(S, K))
      
      P = inter_P / colSums(sweep(theta,MARGIN=c(1),bg, `*`))
      
      F_gradient = array(1, dim = c(N, K))
      alpha = 1
      
      while (sum((alpha * F_gradient)^2) > 0.01) {

        F_gradient = (theta - conv) * sign(conv)
        F_gradient = Multi_F %*% F_gradient
        conv_change = T %*% F_gradient
        temp = T %*% feat / conv_change
        
        alpha = -(max(temp[temp<0]))
        alpha = min(1, alpha * 1.0001)

        feat = feat + (alpha * F_gradient)
        conv = conv_create(N, K, numbase, feat)

      }
      
      LOSS = sum(sweep(P,MARGIN=c(2),colSums(sweep(theta,MARGIN=c(1),bg, `*`)), `*`))
      
      for (i in 1:2) {
        
        inter_theta <- theta[unlist(type[[mid]][i]),]
        theta_dim <- dim(inter_theta)
        inter_theta <- array(inter_theta, dim = c(theta_dim[1],1,K))
        
        inter_P <- array(P, dim = c(S,1,1,K))
        inter_P <- four_recycle(inter_P, 2, theta_dim[1])
        
        inter_mat <- array(mat[i,,], dim = c(1,3,K))
        inter_mat <- three_recycle(inter_mat, 1, theta_dim[1])
        
        temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_theta, `*`)
        temp <- four_recycle(temp, 3, 3)
        temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
        
        inter_bg <- array(bg[unlist(type[[mid]][i])],
                               dim = c(1, length(unlist(type[[mid]][i])), 1))
        inter_bg <- three_recycle(inter_bg, 1, S)
        inter_bg <- three_recycle(inter_bg, 3, 3)
        
        temp <- array(log(sweep(four_colsum(temp, 4),MARGIN=c(1,2,3),inter_bg, `*`)),
                      dim = c(S, length(unlist(type[[mid]][i])), 3))
        
        temp[is.infinite(temp)] <- NA
        LOSS = LOSS - sum(temp * X[,unlist(type[[mid]][i]),], na.rm = TRUE)
        
      }
      
      LOSS = LOSS + 0.5 * reg * sum((theta - conv)^2)
      
      old_LOSS = new_LOSS
      new_LOSS = LOSS
      
    }
    
  }
  
  return("Done")
  
}


