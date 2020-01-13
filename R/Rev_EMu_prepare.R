
# importFrom Biostrings readDNAStringSet countPattern
# mut_count_fast <- function(assembly, mut_file, numbase, nb_chrom) {
#   
#   if(!is.null(assembly) && !is.character(assembly)) {
#     stop("Assembly supplied is not a string.")
#   }
#   
#   genome <- readDNAStringSet(assembly, format = "fasta")
#   
#   if (!is.data.table(mut_file)) {
#     if (!is.data.frame(mut_file) && !is.matrix(mut_file)) {
#       stop("Your input is neither a data.frame or a matrix")
#     }
#     
#     mut_file <- as.data.table(mut_file)
#     
#     mut_file <- mut_file[, chromosome := as.numeric(as.character(chromosome))]
#     mut_file <- mut_file[, chromosome_start := as.numeric(as.character(chromosome_start))]
#     mut_file <- mut_file[, chromosome_end := as.numeric(as.character(chromosome_end))]
#     mut_file <- mut_file[, mutated_from_allele := as.character(mutated_from_allele)]
#     mut_file <- mut_file[, mutated_to_allele := as.character(mutated_to_allele)]
#   }
#   
#   if (numbase != 3 && numbase != 5) {
#     stop("Your specified fragment size is invalid. It should either be 3 or 5")
#   }
#   
#   if (nb_chrom < 1 || nb_chrom > 22) {
#     stop("You specified too little or too many chromosomes for processing")
#   }
#   
#   # Handle mut_matrix sample numbers
#   sample_ids <- unique(mut_file$icgc_sample_id)
#   
#   # How many fragments?
#   if (numbase == 3) {
#     nb_frag = 96
#     background <- rep(0, 96)
#     patterns <- pattern_gen(96)
#     mut_mat <- matrix(0, nrow = length(sample_ids), ncol = 96)
#   } else if (numbase == 5) {
#     nb_frag = 1536
#     background <- rep(0, 1536)
#     patterns <- pattern_gen(1536)
#     mut_mat <- matrix(0, nrow = length(sample_ids), ncol = 1536)
#   }
#   
#   mut_file <- mut_process35(mut_file)
#   
#   # if (fasta == FALSE) {
#   #   
#   #   for (i in 1:nb_chrom) {
#   #     
#   #     chrom_name <- paste("chr", i, sep = "")
#   #     
#   #     for (j in 1:length(patterns)) {
#   #       
#   #       if (length(patterns) == 64) {
#   #         
#   #         temp_inds <- Rfeat2table3(patterns[j]) 
#   #         
#   #       } else if (length(patterns) == 1024) {
#   #         
#   #         temp_inds <- Rfeat2table5(patterns[j]) 
#   #         
#   #       }
#   #       
#   #       found_nb <- countPattern(patterns[j], genome[[chrom_name]])
#   #       
#   #       background[temp_inds[1]] = background[temp_inds[1]] + found_nb
#   #       background[temp_inds[2]] = background[temp_inds[2]] + found_nb
#   #       background[temp_inds[3]] = background[temp_inds[3]] + found_nb
#   #       
#   #     }
#   #     
#   #   }
#   #   
#   # } else if (fasta == TRUE) {
#     
#     chrom_start = 1
#     
#     for (i in 1:nb_chrom) {
#       current_chrom <- genome[[i]]
#       
#       cat("Currently Working on Chromosome ")
#       cat(i)
#       cat("\n")
#       cat("Updating the Genome Background")
#       cat("\n")
#       
#       for (j in 1:length(patterns)) {
#         
#         if (length(patterns) == 64) {
#           
#           temp_inds <- Rfeat2table3(patterns[j]) 
#           
#         } else if (length(patterns) == 1024) {
#           
#           temp_inds <- Rfeat2table5(patterns[j]) 
#           
#         }
#         
#         found_nb <- countPattern(patterns[j], current_chrom)
#         
#         background[temp_inds[1]] = background[temp_inds[1]] + found_nb
#         background[temp_inds[2]] = background[temp_inds[2]] + found_nb
#         background[temp_inds[3]] = background[temp_inds[3]] + found_nb
#         
#       }
#       
#       cat("Updating the Mutation Matrix")
#       cat("\n")
#       
#       if (chrom_start < dim(mut_file)[1]) {
#         # cat("mut_file dimension is:")
#         # cat("\n")
#         # cat(dim(mut_file)[1])
#         # cat("\n")
#         for (k in chrom_start:dim(mut_file)[1]) {
#           # cat("k index is")
#           # cat("\n")
#           # cat(k)
#           # cat("\n")
#           if (mut_file[k]$chromosome > i) {
#             chrom_start = k
#             break
#           } else {
#             
#             if (length(patterns) == 64) {
#               
#               ref_base <- as.character(current_chrom[mut_file[k]$chromosome_start])
#               flankleft_base <- as.character(current_chrom[mut_file[k]$chromosome_start - 1])
#               flankright_base <- as.character(current_chrom[mut_file[k]$chromosome_start + 1])
#               
#               pasted_wind <- paste(flankleft_base, ref_base, flankright_base, sep = "")
#               
#               temp_thisinds <- Rfeat2table3(pasted_wind) 
#               
#               rm(flankleft_base)
#               rm(flankright_base)
#               
#             } else if (length(patterns) == 1024) {
#               
#               ref_base <- as.character(current_chrom[mut_file[k]$chromosome_start])
#               flankleft_base <- as.character(current_chrom[mut_file[k]$chromosome_start - 1])
#               flankleftleft_base <- as.character(current_chrom[mut_file[k]$chromosome_start - 2])
#               flankright_base <- as.character(current_chrom[mut_file[k]$chromosome_start + 1])
#               flankrightright_base <- as.character(current_chrom[mut_file[k]$chromosome_start + 2])
#               
#               pasted_wind <- paste(flankleftleft_base, flankleft_base,
#                                    ref_base, flankright_base, flankrightright_base, sep = "")
#               
#               temp_thisinds <- Rfeat2table5(pasted_wind) 
#               
#               rm(flankleft_base)
#               rm(flankright_base)
#               rm(flankleftleft_base)
#               rm(flankrightright_base)
#               
#             }
#             
#             mut_ref_base <- mut_file[k]$mutated_from_allele
#             
#             trueornot <- ref_base != mut_ref_base
#             
#             if (trueornot) {
#               stop("It seems like the reference genome you supplied is inappropriate")
#             }
#             
#             rm(ref_base)
#             rm(mut_ref_base)
#             rm(trueornot)
#             rm(pasted_wind)
#             
#             # cat("mut_mat dimension is:")
#             # cat("\n")
#             # cat(dim(mut_mat))
#             # cat("\n")
#             
#             incr_loc <- which(sample_ids == mut_file[k]$icgc_sample_id)
#             
#             mut_mat[incr_loc, temp_thisinds[mut_file[k]$mutated_to_allele + 1]] = mut_mat[incr_loc, temp_thisinds[mut_file[k]$mutated_to_allele + 1]] + 1
#             
#             # cat("incr_loc is")
#             # cat("\n")
#             # cat(incr_loc)
#             # cat("\n")
#             # cat("column accessor ir")
#             # cat("\n")
#             # cat(temp_thisinds[mut_file[k]$mutated_to_allele + 1])
#             # cat("\n")
#             
#             rm(incr_loc)
#             
#             chrom_start = k
#           }
#           gc()
#         }
#         
#       }
#     
#       
#     }
#     
#   #}
#   
#   shallowres <- new("Shallowres", mut_mat = mut_mat, wt = background)
#   
#   invisible(shallowres)
#   
# }

Rfeat2table3 <- function(basestring) {
  
  bases <- unlist(strsplit(basestring, ""))
  
  n1 = Base2Num(bases[1])
  n2 = Base2Num(bases[2])
  n3 = Base2Num(bases[3])
  n1t = n1
  
  if (n2 < 2) {
    n1 = 3 - n3
    n3 = 3 - n1t
  } else {
    n2 = 3 - n2
  }
  
  value = (n1*8*3)+(n2*4*3)+(n3*3) + 1
  values <- rep(0, 3)
  values[1] = value
  values[2] = value + 1
  values[3] = value + 2
  
  invisible(values)
  
}

Rfeat2table5 <- function(basestring) {
  
  bases <- unlist(strsplit(basestring, ""))
  
  n1 = Base2Num(bases[1])
  n2 = Base2Num(bases[2])
  n3 = Base2Num(bases[3])
  n4 = Base2Num(bases[4])
  n5 = Base2Num(bases[5])
  
  n1t = n1
  n2t = n2
  
  if (n3 < 2) {
    n1 = 3 - n5
    n2 = 3 - n4
    n4 = 3 - n2t
    n5 = 3 - n1t
  } else {
    n3 = 3 - n3
  }
  
  value = (n1*128*3)+(n2*32*3)+(n3*16*3)+(n4*4*3)+(n5*3) + 1
  values <- rep(0, 3) 
  values[1] = value
  values[2] = value + 1
  values[3] = value + 2
  
  invisible(values)
  
}

pattern_gen <- function(nb_frag) {
  
  all_bases <- c("A", "T", "C", "G")
  
  if (nb_frag == 96) {
    
    patterns <- rep("", 32)
    
    count_keep = 1
    
    for (i in 1:4) {
      for (j in 1:4) {
        for (k in 1:4) {
            
            now_string <- paste(all_bases[i], all_bases[j], all_bases[k], sep = "")
            
            patterns[count_keep] = now_string
            
            count_keep = count_keep + 1
            
        }
      }
    }
    
  } else if (nb_frag == 1536) {
    
    patterns <- rep("", 512)
    
    count_keep = 1
    
    for (i in 1:4) {
      for (j in 1:4) {
        for (k in 1:4) {
          for (q in 1:4) {
            for (p in 1:4) {
              
              now_string <- paste(all_bases[i], all_bases[j], all_bases[k], all_bases[q],
                                  all_bases[p], sep = "")
              
              patterns[count_keep] = now_string
              
              count_keep = count_keep + 1
              
            }
          }
        }
      }
    }
    
  }
  
  invisible(patterns)
  
}

