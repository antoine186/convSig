#' Reads a vcf file and processes it for usage in EM transformations
#'
#' @param datapath The path leading to your multi-sampled vcf file.
#' 
#' @param geno An ID Format string in your vcf file indicating the index of the 
#' genotype information.
#' 
#' @param assembly The path leading to your assembly file (.fa or .fa.gz).
#' 
#' @param five A boolean variable. A value of \code{TRUE} will lead to the function
#' scanning the input files for 5 bases mutation signatures as opposed to 3 bases
#' signatures. a value of \code{FALSE} causes the function to scan for 3 bases signatures.
#'
#' @return A background mutation signatures vector (\code{wt}), which provides
#' the frequency of each possible signature given an assembly file. A matrix (\code{mut_mat})
#' containing the mutational rate of each signature for each sample in your supplied mutation
#' input file.
#' 
#' @examples
#' assembly <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
#' mut_file <- "multisampled_mutations.vcf"
#' 
#' mut_sign <- vcf2mut(mut_file, "GT", assembly)
#' 
#' @export
#' 
#' @importFrom data.table setnames setcolorder
#' @importFrom purrr map2
vcf2mut <- function(datapath, geno = "GT", assembly, five = FALSE) {
  
  # Reading in the vcf file and keep only the SNP
  vcf_data <- read_vcf(datapath)
  
  vcf_data[, `:=` (name = NULL, Vartype = NULL, Qual = NULL, Filter = NULL, Info = NULL)]
  
  # Genotype discovery
  cat("Collecting genotype data from entries in the vcf file")
  cat("\n")
  simp_format <- map2(vcf_data[, FORMAT], geno, format_comprehend)
  vcf_data[, FORMAT := simp_format]
  rm(simp_format)
  
  cnames = colnames(vcf_data)
  minus_start <- grep("FORMAT", cnames)
  end_point <- dim(vcf_data)[2]
  
  cat("Counting the number of mutations for each sample")
  cat("\n")
  cat("This could take a long while")
  cat("\n")
  
  # Hi Anas, this part consumes a lot of memory
  # It takes the FORMAT column (which tells us where GT is) and then takes each sample
  # one column at a time to work on
  # Reminder: I recall us having a discussion about how map2 is even faster than
  # C++
  for (i in c((minus_start+1):end_point)) {
    
    #cur_colname = cnames[i]
    vcf_data[, cnames[i] := map2(vcf_data[, FORMAT], vcf_data[[i]], format_user)] 
    
  }
  
  cat("Mutation count completed")
  cat("\n")
  cat("Formatting intermediate structure")
  cat("\n")
  
  sum_names <- cnames[c((minus_start+1):end_point)]
  chopped_samples <- vcf_data[, sum_names, with = FALSE]
  icgc_height <- sum(unlist(chopped_samples))
  chopped_samples <- data.matrix(chopped_samples)
  
  raw_icgc_form <- icgc_creater(vcf_data[, c(1,2,4,5), with = FALSE], chopped_samples,
                                sum_names, icgc_height, dim(vcf_data)[1])
  
  cat("Filtering and cleaning intermediate structure")
  cat("\n")
  
  proc_icgc_form <- icgc2mut(raw_icgc_form, using.vcf = TRUE)
  
  proc_icgc_form <- icgc_curate(proc_icgc_form)
  
  cat("Mutation signature frequency counting")
  cat("\n")
  
  res <- mut_count(assembly, proc_icgc_form, five = five)
  
  return(res)
  
}

format_comprehend <- function(format_str, geno) {
  
  which_ind = 0
  split_ar <- unlist(strsplit(format_str, ":"))
  
  for (i in c(1:length(split_ar))) {
    
    if (split_ar[i] == geno) {
      
      which_ind = i
      break
    }
    
  }
  
  return(which_ind)
  
}

format_user <- function(pos, sample_inf) {
  
  split_ar <- unlist(strsplit(sample_inf, ":"))
  
  if (grepl("1", split_ar[pos])) {
    
    return(1)
    
  } else {
    
    return(0)
    
  }
  
}



