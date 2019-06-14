#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

complexconv_create <- function(N, K, numbase, feat, mid) {
  
  conv = matrix(1, nrow = N, ncol = K)
  
  for (i in 1:N) {
    temp = i
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

exp_operation <- function(X, bg, conv, P, mat, N, S, K, 
                          type, numbase, bg_test, X_test, mid) {
  
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
    
    
    
  }
  
  return("Done")
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
                         type, numbase, bg_test, X_test, mid)
  
  invisible(exp_res)
  
}
