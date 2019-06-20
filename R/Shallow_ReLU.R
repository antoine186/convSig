#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

#' Performs the ReLU transform on mutational count data
#' 
#' @param mut_obj An object of class 'Shallowres' as produced by the function 
#' \link[=mut_count]{mut_count()}
#' 
#' @param five A boolean or a variable referring to a boolean, which specifies
#'  whether you want to apply a transformation based on a 5-nucleotide 
#'  convolution window
#' 
#' @param K An integer or a variable specifying an integer. This indicates the number
#' of mutational processes you want to detect in your mutational count data via
#'  the transformation.
#'  
#' @return A Feature matrix (\code{feat}), which contains the convolution weights 
#' associated with each mutational processes. An M matrix (\code{mat}), which 
#' contains the probability for all mutation types. A P matrix (\code{P}), which 
#' contains the mutational intensity/activity of each mutational process. A LOSS
#' variable (\code{LOSS}), which displays the LOSS value achieved by the ReLU optimisation.
#' A testing LOSS variable (\code{test_LOSS}), which displays the LOSS value achieved
#' on the testing samples.
#' 
#' @examples
#' relu_res <- relu_transform(EMu_prepped, five = TRUE, K = 6)
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
      inter_bg <- three_recycle(inter_bg, 1, S)
      
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
            
            #inter <- three_colsum((beta_ar / sweep(alpha_array,
                                                   #MARGIN=c(2),x_new, `+`)), 3)
            
            inter <- array(rowSums((beta_ar / 
                                      sweep(alpha_array,MARGIN=c(2),x_new, `+`)),
                                   dims = 2), dim = c(S,N))
            
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
          
          #inter <- three_colsum((beta_ar / sweep(alpha_array,
                                                 #MARGIN=c(2),x_new, `+`)), 3)
          
          inter <- array(rowSums((beta_ar / 
                                    sweep(alpha_array,MARGIN=c(2),x_new, `+`)),
                                 dims = 2), dim = c(S,N))
          
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
        
        if (is.finite(max(temp[temp<0]))) {
        alpha = -(max(temp[temp<0]))
        } else {
          alpha = 1
        }
        
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
  
  test_LOSS = sum(sweep(P,MARGIN=c(2),colSums(sweep(theta,MARGIN=c(1),
                                                    bg_test, `*`)), `*` ))
  
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
    
    inter_bg <- array(bg_test[unlist(type[[mid]][i])],
                      dim = c(1, length(unlist(type[[mid]][i])), 1))
    inter_bg <- three_recycle(inter_bg, 1, S)
    inter_bg <- three_recycle(inter_bg, 3, 3)
    
    temp <- array(log(sweep(four_colsum(temp, 4),MARGIN=c(1,2,3),inter_bg, `*`)),
                  dim = c(S, length(unlist(type[[mid]][i])), 3))
    
    temp[is.infinite(temp)] <- NA
    test_LOSS = test_LOSS - sum(temp * X_test[,unlist(type[[mid]][i]),], na.rm = TRUE)
    
  }
  
  relu_o <- new("relu_obj", feat = feat, mat = mat, P = P,
                LOSS = LOSS, test_LOSS = test_LOSS)
  return(relu_o)
  
}


