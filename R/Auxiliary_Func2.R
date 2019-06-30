reluexp_conv_all <- function(bg, mat, N, K, type, numbase, feat, mid, relu = TRUE) {
  
  if (relu == TRUE) {
    conv <- conv_create(N, K, numbase, feat)
  } else {
    conv <- conv_reverse(N, K, numbase, feat, mid)
  }
  
  conv_all <- array(0, dim = c(N, 3, K))
  
  list_len <- length(unlist(type[[mid]][1]))
  
  for (i in 1:2) {
    
    inter_conv = array(conv[unlist(type[[mid]][i]),], dim = c(list_len, 1, K))
    
    inter_bg = array(bg[unlist(type[[mid]][i])], dim = c(list_len, 1, 1))
    inter_bg <- three_recycle(inter_bg, 3, K)
    
    inter_res <- inter_conv * (inter_bg / sum(bg))
    
    inter_mat = array(mat[i,,], dim = c(1, 3, K))
    inter_mat <- three_recycle(inter_mat, 1, list_len) 
    
    conv_all[unlist(type[[mid]][i]),,] = three_recycle(inter_res, 2, 3) * inter_mat
    
  }
  
  #conv_all <- c(t(as.matrix(conv_all)))
  new_conv_all <- array(0, dim = c(N*3, K))
  
  conv_count = 1
  for (i in 1:N) {
    for (j in 1:3) {
      new_conv_all[conv_count, ] = conv_all[i,j,]
      conv_count = conv_count + 1
    }
  }
  
  invisible(new_conv_all)
  
}

conv_reverse <- function(N, K, numbase, feat, mid) {
  
  conv <- matrix(1, nrow = N, ncol = K)
  
  for (i in 1:N) {
    temp = i - 1
    
    for (j in 1:numbase) {
      if (j == mid) {
        
        conv[i,] = conv[i,] * feat[j,((temp %% 2) + 1),]
        temp = floor(temp/2)
        
        #conv[i,:]*=Feature[j,temp%2,:];
        #temp/=2;
      } else {
        
        conv[i,] = conv[i,] * feat[(numbase - (j - 1) - 1 + 1),((temp %% 4) + 1),]
        temp = floor(temp/4)
        
      }
    }
    
  }
  
  invisible(conv)
  
}
