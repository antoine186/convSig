#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

setClass (
  # Class name
  "exp_obj",
  
  # Defining slot type
  representation (
    feat = "array",
    mat = "array",
    P = "matrix",
    LOSS = "numeric",
    test_LOSS = "numeric"
  )
)

#' Constructs the sparse approximation of h(T_i * F_i) <- check this is true
complexconv_create <- function(N, K, numbase, feat, mid) {
  
  conv = matrix(1, nrow = N, ncol = K)
  
  for (i in 1:N) {
    temp = i - 1
    for (j in 1:numbase) {
      
      if (j == mid) {
        conv[i,] = conv[i,] * feat[j,(temp%%2) + 1,]
        temp = floor(temp/2)
      } else {
        conv[i,] = conv[i,] * feat[j,(temp%%4) + 1,]
        temp = floor(temp/4)
      }
      
    }
  }
  
  return(conv)
}

# A function which helps indexing a 2D matrix using complex list indexing
complexlist_indexer <- function(ar, type, ax, N, K, mid_w8 = FALSE) {
  
  list1 <- unlist(type[[ax]][1])
  list2 <- unlist(type[[ax]][2])
  list3 <- unlist(type[[ax]][3])
  list4 <- unlist(type[[ax]][4])
  
  if (mid_w8 == FALSE) {
    mut_nb <- N/4
    new_ar <- array(0, dim = c(4, mut_nb, K))
  } else {
    mut_nb <- N/2
    new_ar <- array(0, dim = c(2, mut_nb, K))
  }
  
  if (mid_w8 == FALSE) {
    for (i in 1:4) {
      list_mut <- unlist(type[[ax]][i])
      for (j in 1:mut_nb) {
        fancy_ind = list_mut[j]
        for (k in 1:K) {
          new_ar[i,j,k] = ar[fancy_ind, k]
        }
      }
    }
  } else {
    for (i in 1:2) {
      list_mut <- unlist(type[[ax]][i])
      for (j in 1:mut_nb) {
        fancy_ind = list_mut[j]
        for (k in 1:K) {
          new_ar[i,j,k] = ar[fancy_ind, k]
        }
      }
    }
  }
  
  invisible(new_ar)
}

exp_operation <- function(X, bg, conv, P, mat, N, S, K, 
                          type, numbase, bg_test, X_test, mid, Z, feat) {
  
  LOSS = sum(sweep(P,MARGIN=c(2),colSums(sweep(conv,MARGIN=c(1),bg, `*`)), `*`))
  
  for (i in 1:2) {
    
    inter_conv <- conv[unlist(type[[mid]][i]),]
    conv_dim <- dim(inter_conv)
    inter_conv <- array(inter_conv, dim = c(conv_dim[1],1,K))
    
    inter_P <- array(P, dim = c(S,1,1,K))
    inter_P <- four_recycle(inter_P, 2, conv_dim[1])
    
    inter_mat <- array(mat[i,,], dim = c(1,3,K))
    inter_mat <- three_recycle(inter_mat, 1, conv_dim[1])
    
    temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_conv, `*`)
    temp <- four_recycle(temp, 3, 3)
    temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
    temp <- four_colsum(temp, 4)
    
    inter_bg <- array(bg[unlist(type[[mid]][i])], dim = c(1,conv_dim[1],1))
    inter_bg <- three_recycle(inter_bg, 3, 3)
    inter_bg <- three_recycle(inter_bg, 1, S)
    
    temp <- log(sweep(temp,MARGIN=c(1,2,3),inter_bg, `*`))
    
    LOSS = LOSS - sum(sweep(temp,MARGIN=c(1,2,3),X[,unlist(type[[mid]][i]),], `*`))

  }
  
  complete_loss = sum(sweep(P,MARGIN=c(2),colSums(sweep(conv,MARGIN=c(1),bg, `*`)), `*`))
  complete_loss = complete_loss + sum(array(four_colsum((Z * log(Z)), 4),
                                            dim = c(S,N,3)) * X)
  
  for (i in 1:2) {
    
    inter_conv <- conv[unlist(type[[mid]][i]),]
    conv_dim <- dim(inter_conv)
    inter_conv <- array(inter_conv, dim = c(conv_dim[1],1,K))
    
    inter_P <- array(P, dim = c(S,1,1,K))
    inter_P <- four_recycle(inter_P, 2, conv_dim[1])
    
    inter_mat <- array(mat[i,,], dim = c(1,3,K))
    inter_mat <- three_recycle(inter_mat, 1, conv_dim[1])
    
    temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_conv, `*`)
    temp <- four_recycle(temp, 3, 3)
    temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
    
    complete_loss = complete_loss - sum(array(four_colsum(Z[,unlist(type[[mid]][i]),,]
    * log(four_recycle(temp, 4, K)), 4),
    dim = c(S, length(unlist(type[[mid]][i])), 3)) * X[,unlist(type[[mid]][i]),])
    
  }
  
  old_LOSS = 0
  new_LOSS = LOSS
  
  test_LOSS = sum(sweep(P,MARGIN=c(2),
                    colSums(sweep(conv,MARGIN=c(1),bg_test, `*`)), `*`))
  for (i in 1:2) {
    
    inter_conv <- conv[unlist(type[[mid]][i]),]
    conv_dim <- dim(inter_conv)
    inter_conv <- array(inter_conv, dim = c(conv_dim[1],1,K))
    
    inter_P <- array(P, dim = c(S,1,1,K))
    inter_P <- four_recycle(inter_P, 2, conv_dim[1])
    
    inter_mat <- array(mat[i,,], dim = c(1,3,K))
    inter_mat <- three_recycle(inter_mat, 1, conv_dim[1])
    
    temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_conv, `*`)
    temp <- four_recycle(temp, 3, 3)
    temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
    temp_ar <- three_recycle(array(bg_test[unlist(type[[mid]][i])],
                     dim = c(1,length(unlist(type[[mid]][i])),1)), 1, S)
    temp <- log(sweep(four_colsum(temp, 4),MARGIN=c(1,2,3),
                three_recycle(temp_ar, 3, 3), `*`))
    
    test_LOSS = test_LOSS - sum(sweep(temp,MARGIN=c(1,2,3),
                                      X_test[,unlist(type[[mid]][i]),], `*`))
    
  }
  
  while (abs(new_LOSS - old_LOSS) > 10^(1)) {
    
    for (i in 1:2) {
      
      inter_conv <- conv[unlist(type[[mid]][i]),]
      conv_dim <- dim(inter_conv)
      inter_conv <- array(inter_conv, dim = c(conv_dim[1],1,K))
      
      inter_P <- array(P, dim = c(S,1,1,K))
      inter_P <- four_recycle(inter_P, 2, conv_dim[1])
      
      inter_mat <- array(mat[i,,], dim = c(1,3,K))
      inter_mat <- three_recycle(inter_mat, 1, conv_dim[1])
      
      temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_conv, `*`)
      temp <- four_recycle(temp, 3, 3)
      temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
      
      Z[,unlist(type[[mid]][i]),,] = temp
      
    }
    
    Z = Z / array(four_colsum(Z, 4), dim = c(S, N, 3, K))
    
    list_indexer(X, type, mid, N, K)
    list_indexer(Z, type, mid, N, K, X = FALSE)
    
    inter_mat = sweep(list_indexer(Z, type, mid, N, K, X = FALSE),
                      MARGIN=c(1,2,3,4),
                      array(list_indexer(X, type, mid, N, K),
                            dim = c(S,2,length(unlist(type[[mid]][1])),3,1)), `*`)
    mat = array(colSums(five_colsum(inter_mat, 3)), dim = c(2,3,K))
    mat = sweep(mat,MARGIN=c(1,3),three_colsum(mat, 2), `/`)
    
    product = array(colSums(four_colsum(sweep(Z,MARGIN=c(1,2,3),
                                              array(X, dim = c(S,N,3,1)), `*`), 3)),
                    dim = c(N,K))

    for (i in 1:numbase) {

      if (i == mid) {
        feat[i,1:2,] = 1
      } else {
        feat[i,,] = 1
      }

      conv <- complexconv_create(N, K, numbase, feat, mid)

      w8 <- sweep(sweep(conv,MARGIN=c(1),bg, `*`),MARGIN=c(2),colSums(P), `*`)

      if (i != mid) {
        w8_i <- three_colsum(complexlist_indexer(w8, type, i, N, K), 2)
        w8_i <- array(w8_i, dim = c(4, K))
        product_i = three_colsum(complexlist_indexer(product, type, i, N, K), 2)
        product_i <- array(product_i, dim = c(4, K))
      } else {
        w8_i <- three_colsum(complexlist_indexer(w8, type, i,
                                                 N, K, mid_w8 = TRUE), 2)
        w8_i <- array(w8_i, dim = c(2, K))
        product_i = three_colsum(complexlist_indexer(product, type,
                                                     i, N, K, mid_w8 = TRUE), 2)
        product_i <- array(product_i, dim = c(2, K))
      }

      if (i == mid) {
        feat[i,1:2,] = product_i / w8_i
      } else {
        feat[i,,] = product_i / w8_i
      }

      feat[i,,] = sweep(feat[i,,],MARGIN=c(2),
                        array(colSums(feat[i,,]), dim = c(1,K)), `/`)

    }
    
    conv = complexconv_create(N, K, numbase, feat, mid)
    
    step1 <- array(four_colsum(
      four_colsum(sweep(Z,MARGIN=c(1,2,3),X, `*`), 3), 2), dim = c(S,K))
    step2 <- colSums(sweep(conv,MARGIN=c(1),bg, `*`))
    P <- sweep(step1,MARGIN=c(2),step2, `/`) 
    
    old_LOSS = LOSS
    
    LOSS = sum(sweep(P,MARGIN=c(2),colSums(sweep(conv,MARGIN=c(1),bg, `*`)), `*`))
    for (i in 1:2) {
      
      inter_conv <- conv[unlist(type[[mid]][i]),]
      conv_dim <- dim(inter_conv)
      inter_conv <- array(inter_conv, dim = c(conv_dim[1],1,K))
      
      inter_P <- array(P, dim = c(S,1,1,K))
      inter_P <- four_recycle(inter_P, 2, conv_dim[1])
      
      inter_mat <- array(mat[i,,], dim = c(1,3,K))
      inter_mat <- three_recycle(inter_mat, 1, conv_dim[1])
      
      temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_conv, `*`)
      temp <- four_recycle(temp, 3, 3)
      temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
      temp_ar <- three_recycle(array(bg[unlist(type[[mid]][i])],
                                     dim = c(1,length(unlist(type[[mid]][i])),1)), 1, S)
      temp <- log(sweep(four_colsum(temp, 4),MARGIN=c(1,2,3),
                        three_recycle(temp_ar, 3, 3), `*`))
      temp[is.infinite(temp)] <- NA
      
      LOSS = LOSS - sum(sweep(temp,MARGIN=c(1,2,3),
                                        X[,unlist(type[[mid]][i]),], `*`), na.rm = TRUE)
      
    }
    
    new_LOSS = LOSS
    
    cat(new_LOSS - old_LOSS)
    cat("\n")
    
    if (new_LOSS > old_LOSS) {
      stop("Something went wrong during optimization. Probably not something you can solve")
    }
    
    test_LOSS = sum(sweep(P,MARGIN=c(2),
                          colSums(sweep(conv,MARGIN=c(1),bg_test, `*`)), `*`))
    for (i in 1:2) {
      
      inter_conv <- conv[unlist(type[[mid]][i]),]
      conv_dim <- dim(inter_conv)
      inter_conv <- array(inter_conv, dim = c(conv_dim[1],1,K))
      
      inter_P <- array(P, dim = c(S,1,1,K))
      inter_P <- four_recycle(inter_P, 2, conv_dim[1])
      
      inter_mat <- array(mat[i,,], dim = c(1,3,K))
      inter_mat <- three_recycle(inter_mat, 1, conv_dim[1])
      
      temp <- sweep(inter_P,MARGIN=c(2,3,4),inter_conv, `*`)
      temp <- four_recycle(temp, 3, 3)
      temp <- sweep(temp,MARGIN=c(2,3,4),inter_mat, `*`)
      temp_ar <- three_recycle(array(bg_test[unlist(type[[mid]][i])],
                                     dim = c(1,length(unlist(type[[mid]][i])),1)), 1, S)
      temp <- log(sweep(four_colsum(temp, 4),MARGIN=c(1,2,3),
                        three_recycle(temp_ar, 3, 3), `*`))
      
      test_LOSS = test_LOSS - sum(sweep(temp,MARGIN=c(1,2,3),
                                        X_test[,unlist(type[[mid]][i]),], `*`))
      
    }
    
  }
  
  exp_o <- new("exp_obj", feat = feat, mat = mat, P = P, LOSS = LOSS,
               test_LOSS = test_LOSS)
  
  return(exp_o)
}

exp_transform <- function(mut_obj, five = FALSE, K = 5) {
  
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
  
  feat <- array(runif(numbase * 4 * K), dim = c(numbase, 4, K))
  feat[mid,3:4,] = 0
  feat <- sweep(feat,MARGIN=c(1,3),three_colsum(feat, 2), `/`)
  
  conv <- complexconv_create(N, K, numbase, feat, mid)
  
  mat <- array(runif(2*3*K), dim = c(2, 3, K))
  summed_mat <- three_colsum(mat, 2)
  mat <- sweep(mat,MARGIN=c(1,3),summed_mat, `/`)
  
  type <- fragbase_indexer(numbase, N)
  
  exp_res <- exp_operation(X, bg, conv, P, mat, N, S, K, 
                         type, numbase, bg_test, X_test, mid, Z, feat)
  
  invisible(exp_res)
  
}
