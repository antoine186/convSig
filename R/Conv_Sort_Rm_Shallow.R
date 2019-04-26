#' Converts an ICGC file into a Mutation file
#' 
#' @param datapath A string or a variable referencing a string object. This is the path leading to your ICGC file (tsv or csv only). Alternatively, this is your mutation data if you set \code{data.loaded} to \code{TRUE}.
#' @param assembly A string or a variable referencing a string object. This indicates the assembly version used in your genome experiment. Default is set to \code{NULL}, but you really should specify this. \emph{If unspecified, the function will process all of the mutations in your file even if multiple assembly versions are present}.
#' @param Seq A string or a variable referencing a string. This indicates the sequencing strategy/approach used in your genome experiment. Default is set to \code{NULL}, but you really should specify this. \emph{If unspecified, the function will process all of the mutations in your file even if multiple sequencing strategies are present}.
#' @param data.loaded A boolean variable, which indicates whether your input data is already loaded in your environment. The default is set to \code{FALSE} therefore the function looks for a path. If set to \code{TRUE}, the function will work with your data directly in your environment. Make sure that your input data is not malformed; See \code{\link[convSig]{loadICGCexample}} for an example of an acceptable input. Data.frames and matrices are acceptable. Please note that the function is slower with this option.
#'
#' @return A mutation file containing 6 fields/variables: The ICGC sample ID, the chromosome ID, the chromosome start position, the chromosome end position, the reference allele, and the alternate allele.
#' 
#' @section Details:
#' Your input ICGC file must have a header abiding to the ICGC format. The presence of column headers 'mutated_from_allele' and 'mutated_to_allele' are absolute requirements for the usage of this function. You can also technically omit 'assembly_version' and 'sequencing_strategy' if you do not pass them as arguments. However, if you do, then they
#' become requirements. The 'chromosome' header can also be omitted; This is not recommended however as the function will process mutations in the X and Y chromosomes instead of skipping them if they are present in the input file. 
#' 
#' @examples
#' res <- ICGC2Mut("simple_somatic_mutation.open.COCA-CN.tsv")
#' datapath <- "simple_somatic_mutation.open.COCA-CN.tsv"
#' res <- ICGC2Mut(datapath, "GRCh37", "WGS")
#' 
#' @export 
ICGC2Mut <- function(datapath, assembly = NULL, Seq = NULL, data.loaded = FALSE) {
  # Check that data is a string
  if (data.loaded == FALSE) {
    if(!is.character(datapath)) {
      stop("Data file path is not a string.")
    }
  }
  if(!is.null(assembly) && !is.character(assembly)) {
    stop("Assembly supplied is not a string.")
  }
  if(!is.null(Seq) && !is.character(Seq)) {
    stop("Sequencing strategy supplied is not a string.")
  }
  
  # Import ICGC file as a data.table to improve performance
  cat("Loading ICGC file\n")
  tryCatch(
    {
      if (data.loaded == FALSE) {
        x <- data.table::fread(datapath)
      } else {
        x <- data.table::as.data.table(datapath)
      }
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
  # Processing chromosome
  tryCatch(
    {
      x <- x[!chromosome == "X"]
      x <- x[!chromosome == "Y"]
    },
    error=function(cnd) {
      message(paste("Input column 'chromosome' for chromosome ID doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("ICGC2Mut aborted.")
    }
  )
  
  cat("Converting base letters to numbers, this might take a few minutes... \n")
  # Calling Base2Num for mutated_from_allele and mutated_to_allele
  tryCatch(
    {
      from_allele <- purrr::map(x[, mutated_from_allele], Base2Num)
    }, 
    error=function(cnd) {
      message(paste("Input column 'mutated_from_allele' doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("ICGC2Mut aborted.")
    }
  )
  tryCatch(
    {
      to_allele <- purrr::map(x[, mutated_to_allele], Base2Num)
    }, 
    error=function(cnd) {
      message(paste("Input column 'mutated_to_allele' doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("ICGC2Mut aborted.")
    }
  )
  
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
  cat("Returning the result to your specified variable...\n")
  x <- x[from_allele < 4 & to_allele < 4, .(icgc_sample_id, chromosome, 
                                                            chromosome_start, chromosome_end,
                                                            mutated_from_allele, mutated_to_allele)]
  if(data.loaded == TRUE) {
    x <- x[, chromosome := as.numeric(chromosome)]
    x <- x[, chromosome_start := as.numeric(chromosome_start)]
    x <- x[, chromosome_end := as.numeric(chromosome_end)]
    x <- x[, icgc_sample_id := as.character(icgc_sample_id)]
    x <- x[, mutated_from_allele := as.character(mutated_from_allele)]
    x <- x[, mutated_to_allele := as.character(mutated_to_allele)]
  }
  
  cat("Tip: Use data.table::fwrite to write the result to a csv file for example.")
  invisible(x)
}

# Base to number conversion
Base2Num <- function(letter) {
  switch(toupper(letter),
    A = 0,
    G = 1,
    C = 2,
    T = 3,
    4
  )
}

#' Loads an example ICGC file
#' 
#' @return A valid ICGC input file.
#' 
#' @examples
#' example_data <- loadICGCexample()
#' 
#' @export 
loadICGCexample <- function() {
  # Import example input dataset
  load("./Data/example_mutation_dataset.Rda")
  invisible(example_mutation_dataset)
}

# Sorts a mutation file
ICGC_sort <- function(mut_file) {
  cat("Sorting the mutation file\n")
  
  if (all.equal(mut_file[, .(chromosome_start)], mut_file[, .(chromosome_end)], check.attributes = FALSE) != TRUE) {
    warning("Mutations other than single nucleotide mutation detected...\n")
    mut_file <- mut_file[order(chromosome, chromosome_start, chromosome_end)]
  } else {
    mut_file <- mut_file[order(chromosome, chromosome_start)]
  }
  
  invisible(mut_file)
}

# Remove non-single nucleotide polymorphisms
ICGC_snp <- function(mut_file) {
  cat("Removing non-single nucleotide polymorphisms")
  Rcpp::sourceCpp("./src/Shallow_Loops.cpp")
}

# Remove duplicated entries in a mutation file
ICGC_dedup <- function(mut_file) {
  cat("Removing duplicate rows in the mutation file\n")
  mut_file <- unique(mut_file)
  invisible(mut_file)
}


