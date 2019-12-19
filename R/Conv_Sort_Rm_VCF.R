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
vcf2mut <- function(datapath) {
  
  vcf_data <- read_vcf(datapath)
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