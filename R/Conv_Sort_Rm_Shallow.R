#' Converts an ICGC file into a Mutation file
#' 
#' @param datapath A string or a variable referencing a string object. This is the path leading to your ICGC file.
#' @param assembly A string or a variable referencing a string object. This indicates the assembly version used in your genome experiment. Default is set to NULL, but you really should specify this. If unspecified, the function will process all of the mutations in your file even if multiple assembly versions are present.
#' @param Seq A string or a variable referencing a string. This indicates the sequencing strategy/approach used in your genome experiment. Default is set to NULL, but you really should specify this. If unspecified, the function will process all of the mutations in your file even if multiple sequencing strategies are present.
#'
#' @return A mutation file containing 6 fields/variables: The ICGC sample ID, the chromosome ID, the chromosome start position, the chromosome end position, the reference allele, and the alternate allele
#' 
#' @section Details:
#' Your input ICGC file must have a header abiding to the ICGC format. The presence of column headers 'mutated_from_allele' and 'mutated_to_allele' are absolute requirements for the usage of this function. You can also technically omit 'assembly_version' and 'sequencing_strategy' if you do not pass them as arguments. However, if you do, then they
#' become requirements. The 'chromosome' header can also be omitted; This is not recommended however as the function will process mutations in the X and Y chromosomes instead of skipping them if they are present in the input file. 
#' 
#' @examples
#' ICGC2Mut("simple_somatic_mutation.open.COCA-CN.tsv")
#' datapath <- "simple_somatic_mutation.open.COCA-CN.tsv"
#' ICGC2Mut(datapath, "GRCh37", "WGS")
#' 
#' @export 
ICGC2Mut <- function(datapath, assembly = NULL, Seq = NULL) {
  # Check that data is a string
  if(!is.character(datapath)) {
    stop("Data file path is not a string.")
  }
  if(!is.character(assembly)) {
    stop("Assembly supplied is not a string.")
  }
  if(!is.character(Seq)) {
    stop("Sequencing strategy supplied is not a string.")
  }
  
  # Import ICGC file as a data.table to improve performance
  cat("Loading ICGC file\n")
  tryCatch(
    {
      x <- data.table::fread(datapath)
    }, 
    error=function(cnd) {
      message(paste("There seems to be a problem with the filepath you supplied."))
      message("Here's the original error message:")
      message(cnd)
      stop("ICGC2Mut aborted.")
    }
  )
  
  cat("Checking arguments supplied with call\n")
  # Checking for assembly_version
  if (!is.null(assembly)) {
    if (!any(colnames(x) == "assembly_version")) {
      stop("Cannot find a header variable corresponding to assembly_version.\n
          No way to tell which column refers to assembly version.")
    }
    x <- x[assembly_version == assembly]
    if (dim(x)[1] == 0) {
      stop("Specified assembly_version does not exist in the data.")
    }
  } else {
    warning("You should supply a valid assembly version.")
  }
  # Checking for sequencing_strategy
  if (!is.null(Seq)) {
    if (!any(colnames(x) == "sequencing_strategy")) {
      stop("Cannot find a header variable corresponding to sequencing_strategy.\n
          No way to tell which column refers to sequencing strategy.")
    }
    x <- x[sequencing_strategy == Seq]
    if (dim(x)[1] == 0) {
      stop("Specified sequencing_strategy does not exist in the data.")
    }
  } else {
    warning("You should supply a valid sequencing strategy.")
  } 
  # Checking chromosome - Processing chromosome
  if (!any(colnames(x) == "chromosome")) {
    warning("Cannot find a header variable corresponding to chromosome.\n
            We might inadvertently process the X and Y chromosomes!")
  } else {
    x <- x[!chromosome == "X"]
    x <- x[!chromosome == "Y"]
  }
  
  cat("Converting base letters to numbers, this might take a few minutes... \n")
  # Calling Base2Num for mutated_from_allele and mutated_to_allele
  from_allele <- purrr::map(x[, mutated_from_allele], Base2Num)
  to_allele <- purrr::map(x[, mutated_to_allele], Base2Num)
  
  cat("Conversion successful\n
      ´´´´´´¶¶¶¶
      ´´´´¶¶´´´´´¶
      ´´´´´¶´´´´´¶
      ´´´´´´¶´´´´¶
      ´´´´´´¶´´´¶
      ´´´´¶¶¶¶¶¶¶¶¶¶¶¶
      ´´´¶´´´´´´´´´´´´¶
      ´´¶´´´´´´´´´´´´¶
      ´¶¶´´´¶¶¶¶¶¶¶¶¶¶¶
      ´¶´´´´´´´´´´´´´´´¶
      ´¶´´´´´´´´´´´´´´´¶
      ´´¶´´´¶¶¶¶¶¶¶¶¶¶¶
      ´´´¶´´´´´´´´´´´¶
      ´´´´¶¶¶¶¶¶¶¶¶¶¶\n")
  cat("Returning the result...\n")
  x <- x[from_allele < 4 & to_allele < 4, .(icgc_sample_id, chromosome, 
                                                            chromosome_start, chromosome_end,
                                                            mutated_from_allele, mutated_to_allele)]
  cat("Tip: Use data.table::fwrite to write the result to a csv file for example.")
  invisible(x)
}

Base2Num <- function(letter) {
  switch(toupper(letter),
    A = 0,
    G = 1,
    C = 2,
    T = 3,
    4
  )
}


