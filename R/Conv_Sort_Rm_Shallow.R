#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

#' Converts an ICGC file into a mutation file
#' 
#' @param datapath A string or a variable referencing a string object. This is 
#' the path leading to your ICGC file (tsv or csv only). Alternatively, this 
#' is your mutation data if you supplied either a \code{data.frame} or a \code{matrix}.
#' Please make sure that your input data is not malformed; 
#' See \link[=loadICGCexample]{loadICGCexample()} for an example of an 
#' acceptable input. Please note that the function is slower with either a supplied
#' \code{data.frame} or \code{matrix}.
#' @param assembly A string or a variable referencing a string object. This 
#' indicates the assembly version used in your experiment. Default is 
#' set to \code{NULL}, but you really should specify this. \emph{If unspecified,
#' the function will process the most represented assembly version present in your
#' input file}.
#' @param Seq A string or a variable referencing a string. This indicates the 
#' sequencing strategy/approach used in your genome experiment. Default is set 
#' to \code{NULL}, but you really should specify this. \emph{If unspecified, 
#' the function will process the most represented sequencing strategy present in your
#' input file}.
#'
#' @return A mutation file containing 6 fields/variables: The ICGC sample ID, 
#' the chromosome ID, the chromosome start position, the chromosome end 
#' position, the reference allele, and the alternate allele.
#' 
#' @section Details:
#' Your input ICGC file must have a header abiding to the ICGC format. The 
#' presence of column headers 'mutated_from_allele' and 'mutated_to_allele' 
#' are absolute requirements for the usage of this function. You can also 
#' technically omit 'assembly_version' and 'sequencing_strategy' if you do not 
#' pass them as arguments. However, if you do, then they become requirements. 
#' The 'chromosome' header can also be omitted; This is 
#' not recommended however as the function will process mutations in the X and 
#' Y chromosomes instead of skipping them if they are present in the input file. 
#' See \link[=loadICGCexample]{loadICGCexample()} for a closer look at a legal 
#' input format.
#' 
#' @examples
#' res <- icgc2mut("simple_somatic_mutation.open.COCA-CN.tsv")
#' # Or
#' datapath <- "simple_somatic_mutation.open.COCA-CN.tsv"
#' res <- icgc2mut(datapath, "GRCh37", "WGS")
#' 
#' @export 
#' 
#' @importFrom data.table fread data.table as.data.table
#' @importFrom purrr map
icgc2mut <- function(datapath, assembly = NULL, Seq = NULL) {
  
  if(!is.null(assembly) && !is.character(assembly)) {
    stop("Assembly supplied is not a string.")
  }
  if(!is.null(Seq) && !is.character(Seq)) {
    stop("Sequencing strategy supplied is not a string.")
  }
  
  # Set data.loaded accordingly to flag the need for column conversion
  data.loaded = FALSE
  
  # Import ICGC file as a data.table to improve performance
  cat("Loading ICGC file\n")
  tryCatch(
    {
      if ("character" %in% class(datapath)) {
        x <- fread(datapath)
      } else {
        x <- as.data.table(datapath)
        data.loaded = TRUE
      }
    }, 
    error=function(cnd) {
      message(paste("There seems to be a problem with the filepath you supplied."))
      message("Here's the original error message:")
      message(cnd)
      stop("icgc2mut aborted.")
    }
  )
  
  cat("Checking arguments supplied with call\n")
  # Checking for assembly_version
  if (!any(colnames(x) == "assembly_version")) {
    stop("Cannot find a header variable corresponding to assembly_version.\n
          No way to tell which column refers to assembly version.")
  }
  
  if (!is.null(assembly)) {
    x <- x[assembly_version == assembly]
    if (dim(x)[1] == 0) {
      stop("Specified assembly_version does not exist in the data.")
    }
  } else {
    warning("You should supply a valid assembly version.")
    assembly_stat <- x[, .(.N), by = .(assembly_version)]
    assembly_stat <- assembly_stat[order(-N)]
    assembly_stat$assembly_version <- as.character(assembly_stat$assembly_version)
    assembly <- as.character(assembly_stat[1,1])
    
    cat(paste("We have chosen the following assembly:", assembly, "\n", sep = " "))
    cat("\n")
    print(assembly_stat)
    cat("\n")
    
    x <- x[assembly_version == assembly]
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
    seq_stat <- x[, .(.N), by = .(sequencing_strategy)]
    seq_stat <- seq_stat[order(-N)]
    seq_stat$sequencing_strategy <- as.character(seq_stat$sequencing_strategy)
    Seq <- as.character(seq_stat[1,1])
    
    cat(paste("We have chosen the following sequencing strategy:", Seq, "\n", sep = " "))
    cat("\n")
    print(seq_stat)
    cat("\n")
    
    x <- x[sequencing_strategy == Seq]
  } 
  # Processing chromosome
  tryCatch(
    {
      x <- x[!chromosome == "X"]
      x <- x[!chromosome == "Y"]
    },
    error=function(cnd) {
      message(paste("Input column 'chromosome' for chromosome ID doesn't seem 
                    to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("icgc2mut aborted.")
    }
  )
  
  cat("Converting base letters to numbers, this might take a few minutes... \n")
  # Calling Base2Num for mutated_from_allele and mutated_to_allele
  tryCatch(
    {
      from_allele <- map(x[, mutated_from_allele], Base2Num)
    }, 
    error=function(cnd) {
      message(paste("Input column 'mutated_from_allele' doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("icgc2mut aborted.")
    }
  )
  tryCatch(
    {
      to_allele <- map(x[, mutated_to_allele], Base2Num)
    }, 
    error=function(cnd) {
      message(paste("Input column 'mutated_to_allele' doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("icgc2mut aborted.")
    }
  )
  
  cat("Conversion successful\n")
  cat("Returning the result to your specified variable...\n")
  x <- x[from_allele < 4 & to_allele < 4, .(icgc_sample_id, chromosome, 
                                                            chromosome_start, 
                                            chromosome_end, mutated_from_allele,
                                            mutated_to_allele)]
  if(data.loaded == TRUE) {
    x <- x[, chromosome := as.numeric(as.character(chromosome))]
    x <- x[, chromosome_start := as.numeric(as.character(chromosome_start))]
    x <- x[, chromosome_end := as.numeric(as.character(chromosome_end))]
    x <- x[, icgc_sample_id := as.character(icgc_sample_id)]
    x <- x[, mutated_from_allele := as.character(mutated_from_allele)]
    x <- x[, mutated_to_allele := as.character(mutated_to_allele)]
  }
  
  cat("Tip: Use data.table::fwrite to write the result to a csv file for example.\n")
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
#' @return A valid ICGC input file (that can for example be used for \link[=icgc2mut]{icgc2mut()}).
#' 
#' @examples
#' loadICGCexample()
#' 
#' @export 
loadICGCexample <- function() {
  # Import example input dataset
  load("./data/example_mutation_dataset.Rda", envir = parent.frame())
}

# Sorts the mutation file
icgc_sort <- function(mut_file) {
  cat("Sorting the mutation file\n")
  
  mut_file <- mut_file[, chromosome := as.numeric(as.character(chromosome))]
  mut_file <- mut_file[, chromosome_start := as.numeric(as.character(chromosome_start))]
  mut_file <- mut_file[, chromosome_end := as.numeric(as.character(chromosome_end))]
  mut_file <- mut_file[, icgc_sample_id := as.character(icgc_sample_id)]
  mut_file <- mut_file[, mutated_from_allele := as.character(mutated_from_allele)]
  mut_file <- mut_file[, mutated_to_allele := as.character(mutated_to_allele)]
  
  tryCatch(
    {
      if (all.equal(mut_file[, .(chromosome_start)], mut_file[, .(chromosome_end)],
                    check.attributes = FALSE) != TRUE) {
        cat("Mutations other than single nucleotide mutation detected...\n")
        mut_file <- mut_file[order(chromosome, chromosome_start, chromosome_end)]
      } else {
        mut_file <- mut_file[order(chromosome, chromosome_start)]
      }
    }, 
    error=function(cnd) {
      message(paste("It seems that your input file does not conform with the 
                    required format."))
      message("Here's the original error message:")
      message(cnd)
      stop("icgc_curate aborted.")
    }
  )
  
  invisible(mut_file)
}

# Remove non-single nucleotide polymorphisms
icgc_snp <- function(mut_file) {
  cat("Removing non-single nucleotide polymorphisms\n")
  tryCatch(
    {
      rm_ind <- RM_nonSNP(mut_file[,.(chromosome_start, chromosome_end)],
                          rep(TRUE, dim(mut_file)[1]))
      mut_file <- mut_file[rm_ind]
    },
    error=function(cnd) {
      message(paste("Oops, it seems like an internal problem happened. 
                    Probably something you cannot fix on your end."))
      message("Here's the original error message:")
      message(cnd)
      stop("Operation aborted.")
    }
  )
  
  if (dim(mut_file)[1] == 0) {
    stop("Your input file only contains non single nucleotide changes.", 
         "This package only deals with those unfortunately")
  }
  
  invisible(mut_file)
}

# Remove duplicated entries in a mutation file
icgc_dedup <- function(mut_file) {
  cat("Removing duplicate rows in the mutation file\n")
  mut_file <- unique(mut_file)
  invisible(mut_file)
}

#' Sorts and removes duplicate entries in a mutation file output 
#' by \code{icgc2mut()}
#' 
#' @param mut_file A string or a variable referencing a mutation file output 
#' by \link[=icgc2mut]{icgc2mut()}. You can specify a file you created, however
#'  you have to make sure that it has the correct format (i.e. please view 
#'  \emph{Details}). Your input should either be a \code{data.frame}, a 
#'  \code{matrix}, or a \code{data.table}.
#' @param remove.nonSNP A boolean variable indicating whether the function 
#' should remove non-single nucleotide changes. This is set by default 
#' to \code{TRUE} as our downstream filtering method only handles single 
#' nucleotide changes. You can toggle this parameter to \code{FALSE} for your 
#' own purposes, although it should remain \code{TRUE} if you want to utilize 
#' our entire pipeline for your analysis.
#'
#' @return A sorted mutation file with duplicate entries removed and, depending 
#' on user specification, non-single nucleotide changes removed.
#' 
#' @section Details:
#' Your input mutation file must at least have the following column headers: 
#' icgc_sample_id (i.e. the ICGC sample ID), chromosome (i.e. the Chromosome ID),
#'  chromosome_start (i.e. the chromosome 
#' start position), chromosome_end (i.e. the chromosome end position), 
#' mutated_from_allele (i.e. the reference allele), and mutated_to_allele 
#' (i.e. alternate allele).
#' 
#' @examples
#' res <- icgc_curate(mutation_data)
#' 
#' @export
#' @importFrom data.table is.data.table as.data.table
icgc_curate <- function(mut_file, remove.nonSNP = TRUE) {
  
  tryCatch(
    {
      if (!is.data.table(mut_file)) {
        if (!is.data.frame(mut_file) && !is.matrix(mut_file)) {
          stop("Your input is neither a data.frame or a matrix")
        }
        
        if (!any(colnames(mut_file) == "chromosome")) {
          stop("Cannot find a header variable corresponding to \'chromosome\'")
        } else if (!any(colnames(mut_file) == "chromosome_start")) {
          stop("Cannot find a header variable corresponding to \'chromosome_start\'")
        } else if (!any(colnames(mut_file) == "chromosome_end")) {
          stop("Cannot find a header variable corresponding to \'chromosome_end\'")
        }
        
        mut_file <- as.data.table(mut_file)
        
        mut_file <- mut_file[, chromosome := as.numeric(as.character(chromosome))]
        mut_file <- mut_file[, chromosome_start := as.numeric(as.character(chromosome_start))]
        mut_file <- mut_file[, chromosome_end := as.numeric(as.character(chromosome_end))]
        mut_file <- mut_file[, mutated_from_allele := as.character(mutated_from_allele)]
        mut_file <- mut_file[, mutated_to_allele := as.character(mutated_to_allele)]
      }
    }, 
    error=function(cnd) {
      message(paste("Something has gone wrong with the input you supplied..."))
      message("Here's the original error message:")
      message(cnd)
      stop("Operation aborted.")
    }
  )
  
  mut_file <- icgc_sort(mut_file)
  
  if (remove.nonSNP == TRUE) {
    mut_file <- icgc_snp(mut_file)
  } 
  
  mut_file <- icgc_dedup(mut_file)
  invisible(mut_file)
}

chrom_check <- function(interdata) {
  if (!any(colnames(interdata) == "chromosome")) {
    stop("Cannot find a column header by the name of \'chromosome\'")
  }
  
  if ("integer" %in% class(interdata$chromosome)) {
    return(interdata)
  } else if ("character" %in% class(interdata$chromosome)) {
    return("lol")
  } else if ("factor" %in% class(interdata$chromosome)) {
    return("lol")
  }
}

chrom_string_split <- function(chrom) {
  return("lol")
}


