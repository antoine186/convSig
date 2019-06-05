#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

setClass (
  # Class name
  "Shallowres",
  
  # Defining slot type
  representation (
    mut_mat = "matrix",
    wt = "numeric"
  )
)

#' Read in the reference genome fasta file
#'
#' @importFrom data.table fread data.table as.data.table
readfast <- function(datapath) {
  tryCatch(
    {
      #if (!file.exists(datapath)) {
        #stop(datapath, ": File doesn't exist.")
      #} 
      if ("character" %in% class(datapath)) {
        x <- fread(datapath, header = FALSE)
      } else if ("data.frame" %in% class(datapath)) {
        x <- as.data.table(datapath)
        #data.loaded = TRUE
      }
    },
    error=function(cnd) {
      message(paste("There seems to be a problem with the filepath you supplied."))
      message("Here's the original error message:")
      message(cnd)
      stop("Reading in assembly file input aborted.")
    }
  )
  invisible(x)
  
  # Handle factors if we are dealing with data in env.
}

#' Read in the mutation input file
#'
#' @importFrom data.table fread data.table as.data.table
#' 
#' @export
readmut <- function(datapath) {
  tryCatch(
    {
      if ("character" %in% class(datapath)) {
        x <- fread(datapath)
      } else {
        x <- as.data.table(datapath)
        #data.loaded = TRUE
        x <- x[, icgc_sample_id := as.character(icgc_sample_id)]
        x <- x[, chromosome := as.numeric(as.character(chromosome))]
        x <- x[, chromosome_start := as.numeric(as.character(chromosome_start))]
        x <- x[, mutated_from_allele := as.character(mutated_from_allele)]
        x <- x[, mutated_to_allele := as.character(mutated_to_allele)]
      }
    }, 
    error=function(cnd) {
      message(paste("There seems to be a problem with the filepath you supplied."))
      message("Here's the original error message:")
      message(cnd)
      stop("Reading in mutation file input aborted.")
    }
  )
  invisible(x)
}

#' Computes the frequency of each possible trinucleotide or 5-nucleotide mutation signature
#' 
#' @param reference A string or a variable referencing a string object. This is
#' the path leading to your assembly file (.fa or .fa.gz). 
#' 
#' @param mut_file A string or a variable referencing a string object. This is 
#' the path leading to your mutation input file (tsv or csv only). Alternatively,
#' this is your mutation data if you supplied either a \code{data.frame} or a \code{matrix}.
#' Please make sure that your input data is not malformed; It should contain 
#' the following columns: icgc_sample_id (i.e. the ICGC sample ID), chromosome 
#' (i.e. the Chromosome ID), chromosome_start (i.e. the chromosome start position), 
#' mutated_from_allele (i.e. the reference allele), and mutated_to_allele 
#' (i.e. alternate allele). 
#' Your input file must also be sorted and must only contain single nucleotide
#' changes (see \link[=icgc_curate]{icgc_curate()}).
#' \emph{Such files can be produced by some of our functions
#' like \link[=icgc2mut]{icgc2mut()}}.
#' 
#' @param five A boolean variable. A value of \code{TRUE} will lead to the function
#' scanning the input files for 5 bases mutation signatures as opposed to 3 bases
#' signatures. a value of \code{FALSE} causes the function to scan for 3 bases signatures. 
#' This value is set to \code{FALSE} by default.
#' 
#' @return A background mutation signatures vector (\code{wt}), which provides
#' the frequency of each possible signature given an assembly file. A matrix (\code{mut_mat}) 
#' containing the mutational rate of each signature for each sample in your supplied mutation
#' input file.
#' 
#' @section Details: Explain the return value a bit more
#' 
#' @examples 
#' assembly <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
#' mut_file <- "mutation_file_input.tsv"
#' 
#' mut_sign <- mut_count3(assembly, mut_file)
#'
#' @export
mut_count <- function(reference, mut_file, five = FALSE) {
  # Make sure mutation input is a data.table

  cat("Loading the assembly\n")
  # Call readfast
  reference_gen <- readfast(reference)
  
  cat("Loading the mutation input file\n")
  # Read in the mutation input file
  datapath <- readmut(mut_file)
  
  cat("Miscellaneous processing...\n")
  nb_uniq = length(unique(as.character(datapath$icgc_sample_id)))
  # This is the row names of init_mut_mat
  uniq_sample <- unique(as.character(datapath$icgc_sample_id))
  
  if (five == FALSE) {
    init_mut_mat <- matrix(0L, nrow = nb_uniq, ncol = 96) 
    init_wt <- rep(0, 96)
  } else if (five == TRUE) {
    init_mut_mat <- matrix(0L, nrow = nb_uniq, ncol = 1536) 
    init_wt <- rep(0, 1536)
  }
  
  # Call mut_process3 here and store result in variable called treated_mut
  cat("Processing and encoding the mutation input file\n")
  # Order of this file's column is super important <= Here
  treated_mut <- mut_process3(datapath)
  
  cat("Counting the frequency of mutation fragment types. This could take a few minutes...\n")
  shallowres <- new("Shallowres", mut_mat = init_mut_mat, wt = init_wt)
  
  if (five == FALSE) {
    shallowres = shallow_loop3(shallowres, reference_gen, treated_mut, uniq_sample)
  } else if (five == TRUE) {
    shallowres = shallow_loop5(shallowres, reference_gen, treated_mut, uniq_sample)
  }
  
  invisible(shallowres)
}

#' A function that treats the mutation input file and reorders its columns
#' 
#' @importFrom data.table data.table as.data.table
#' @importFrom purrr map
mut_process3 <- function(datapath) {
  
  tryCatch(
    {
      datapath <- datapath[!chromosome == "X"]
      datapath <- datapath[!chromosome == "Y"]
    },
    error=function(cnd) {
      message(paste("Input column 'chromosome' for chromosome ID doesn't seem 
                    to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("###### aborted.")
    }
  )
  
  tryCatch(
    {
      from_allele <- map(datapath[, mutated_from_allele], Base2Num)
    }, 
    error=function(cnd) {
      message(paste("Input column 'mutated_from_allele' doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("####### aborted.")
    }
  )
  tryCatch(
    {
      if ("numeric" %in% class(datapath[, mutated_to_allele])) {
        stop("Either your input mutation file has already been processed 
             or the values in your 'mutated_to_allele' column are not letter
             bases (i.e. wrong format)")
      } else {
        to_allele <- map(datapath[, mutated_to_allele], Base2Num)
      }
    }, 
    error=function(cnd) {
      message(paste("There seems to be a problem with your 
                    'mutated_to_allele' input column."))
      message("Here's the original error message:")
      message(cnd)
      stop("######## aborted.")
    }
  )
  
  datapath <- datapath[, mutated_to_allele := as.numeric(to_allele)]

  datapath <- datapath[from_allele < 4 & 
                        mutated_to_allele < 4, .(icgc_sample_id, chromosome, 
                                                  chromosome_start, 
                                                  mutated_from_allele,
                                                  mutated_to_allele)]
  
  both_alleles <- data.table(from_allele = from_allele, to_allele = to_allele)
  both_alleles <- both_alleles[from_allele < 4 & to_allele < 4, .(from_allele,
                                                                  to_allele)]
  
  
  fromto_allele <- alt_reverse(datapath, as.numeric(both_alleles[, from_allele]))
  datapath <- datapath[, mutated_to_allele := as.numeric(fromto_allele$to_allele)]
  
  invisible(datapath)
}

#' @importFrom data.table data.table as.data.table
alt_reverse <- function(datapath, from_allele) {
  alleles = data.table(from_allele = from_allele, to_allele = 
                         datapath[,mutated_to_allele])
  tryCatch(
    {
      alleles <- reverse_transform(alleles)
    },
    error=function(cnd) {
      message(paste("Input column 'mutated_from_allele' doesn't seem to exist."))
      message("Here's the original error message:")
      message(cnd)
      stop("###### aborted.")
    }
  )
  
  invisible(alleles)
}

#chrom_separate <- function() {
  
#}


