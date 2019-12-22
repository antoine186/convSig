#' Read vcf files and processes them / User correct assembly
#'
#' @param datapath
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom data.table setnames setcolorder
#' @importFrom purrr map2
vcf2mut <- function(datapath, geno = "GT") {
  
  # Reading in the vcf file and keep only the SNP
  vcf_data <- read_vcf(datapath)
  
  vcf_data[, `:=` (name = NULL, Vartype = NULL, Qual = NULL, Filter = NULL, Info = NULL)]
  
  # Genotype discovery
  simp_format <- map2(vcf_data[, FORMAT], geno, format_comprehend)
  vcf_data[, FORMAT := simp_format]
  rm(simp_format)
  
  cnames = colnames(vcf_data)
  minus_start <- grep("FORMAT", cnames)
  end_point <- dim(vcf_data)[2]
  
  for (i in c((minus_start+1):end_point)) {
    
    #cur_colname = cnames[i]
    vcf_data[, cnames[i] := map2(vcf_data[, FORMAT], vcf_data[[i]], format_user)] 
    
  }
  
  sum_names <- cnames[c((minus_start+1):end_point)]
  chopped_samples <- unlist(vcf_data[, sum_names, with = FALSE])
  icgc_height <- sum(unlist(vcf_data[, sum_names, with = FALSE]))
  
  raw_icgc_form <- icgc_creater(vcf_data[, c(1,2,4,5), with = FALSE], chopped_samples,
                                sum_names, icgc_height)
  
  ### Under development
  vcf_data <- vcf_data[, c("Chrom", "Pos", "ID", "REF", "ALT"), with=FALSE]
  
  setnames(vcf_data, c("chromosome", "chromosome_start", "ID",
           "mutated_from_allele", "mutated_to_allele"))
  
  dummy_icgc_col <- rep("dummy_sample", dim(vcf_data)[1])
  
  vcf_data <- vcf_data[, ID:=NULL]
  vcf_data <- vcf_data[, ("chromosome_end") := vcf_data[, "chromosome_start", with=FALSE]]
  vcf_data <- vcf_data[, ("icgc_sample_id") := dummy_icgc_col]
  
  setcolorder(vcf_data, c("icgc_sample_id", "chromosome", "chromosome_start",
                          "chromosome_end", "mutated_from_allele", "mutated_to_allele"))
  
  return(vcf_data)
  
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



