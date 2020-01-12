#' Read in the vcf file
#'
#' @section Details:
#' This is an auxiliary function.
#'
#' @importFrom data.table fread
#' @importFrom dplyr if_else %>%
#' @importFrom stringr str_c str_sub str_detect str_replace str_extract
read_vcf <- function(vcf_file) {
  
  # Reading in the vcf file
  variant_dt <- fread(vcf_file, sep = "\t", skip = "CHROM", fill = TRUE)
  
  # Deleting spurrious columns
  cnames = colnames(variant_dt)
  rm_col = character()
  
  for (i in c(1:length(cnames))) {
    
    if (grepl("^V", cnames[i])) {
      
      rm_col = c(rm_col, cnames[i])
      
    }
    
  }
  
  for (i in c(1:length(rm_col))) {
    
    variant_dt[, rm_col[i] := NULL]
    
  }
   
  rename_vec <- c("Chrom", "Pos", "ID", "REF", "ALT", "Qual", "Filter", "Info")
  for (i in c(1:8)) {
    
    colnames(variant_dt)[i] = rename_vec[i]
    
  }
  
  variant_dt[, `:=`(name = str_c(Chrom, Pos, str_replace(REF, ",", "-"),
                                 str_replace(ALT, ",", "-"), sep = "-"),
                    ref_n = str_c(REF, ","),
                    alt_n = str_c(ALT, ","))
             ][, `:=`(ref_n = nchar(ref_n %>% str_extract(".*/?(?=,)")),
                      alt_n = nchar(alt_n %>% str_extract(".*/?(?=,)")),
                      Vartype = "NA")
               ][, Vartype := if_else(ref_n == 1 & alt_n == 1, "snp", Vartype)
                 ][, Vartype := if_else(ref_n == 1 & alt_n > 1, "ins", Vartype)
                   ][, Vartype := if_else(ref_n > 1 & alt_n == 1, "del", Vartype)
                     ][, Vartype := if_else(ref_n > 2 & alt_n > 2, "sub", Vartype)
                       ][, c("ref_n", "alt_n") := NULL]
  
  variant_dt <- variant_dt[Vartype == "snp"]
  
  return(variant_dt)

}