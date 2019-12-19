#' Read vcf files
#'
#' @param vcf_file
#' @param well_id
#' @param region
#' @param drop_col
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom data.table fread
#' @importFrom VariantAnnotation scanVcfHeader
#' @importFrom VariantAnnotation samples
#' @importFrom dplyr if_else %>%
#' @importFrom stringr str_c str_sub str_detect str_replace str_extract
read_vcf <- function(vcf_file, well_id = NULL, region = NULL,
                     drop_col = FALSE) {
  
  # Reading in the vcf file
  variant_dt <- fread(vcf_file, sep = "\t", skip = "CHROM", fill = TRUE)
  
  # Hidden option to drop columns beyond the mandatory columns
  if (drop_col) {
    cat("\nDropping additional columns found in .vcf file\n")
    variant_dt <- variant_dt[, 1:(9 + length(well_id))]
  }
  
  if (is.null(well_id)) {
    well_id <- scanVcfHeader(vcf_file) %>%
      samples()
  } else if (class(well_id) == "numeric") {
    well_id <- scanVcfHeader(vcf_file) %>%
      samples() %>%
      substr(1, well_id)
  } else if ((length(well_id) + 9 != ncol(variant_dt)) && !drop_col) {
    stop("\n\nNumber of columns do not match number of IDs provided.",
         "\nExpecting: ", ncol(variant_dt) - 9, "\t Provided: ",
         length(well_id), "\n")
  }
  
  colnames(variant_dt) <-
    c("Chrom", "Pos", "ID", "REF", "ALT", "Qual", "Filter", "Info")
  
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
  
  # TODO: including regions          xcv
  variant_dt <- variant_dt[, c("Chrom", "Pos", "name", "ID", "REF", "ALT",
                               "Vartype", "Qual", "Filter", "Info"), with = F]
  
  variant_dt <- variant_dt[Vartype == "snp"]
  
  return(variant_dt)

}