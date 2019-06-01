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

#res <- conv_create(32, 5, 3)
#res <- conv_create(512, 5, 5)

# "/Users/antoinetian/Documents/GitHub/convSig-shallow/feat3test.csv"

# P <- matrix(runif(2*2), nrow = 2, ncol = 2)
# P2 <- matrix(runif(2*2), nrow = 2, ncol = 2)
# 
# bg <- array(runif(2*2*2), c(2, 2, 2))
# theta <- array(runif(2*2*3), c(2, 2, 2))
# P <- array(runif(10*5), c(10, 5))
# 
# int <- base::colSums(bg * theta)
# 
# res <- P%*%P2
# lol <- bg%*%theta



