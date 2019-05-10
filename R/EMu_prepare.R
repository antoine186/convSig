#' @useDynLib convSig
#' @importFrom Rcpp sourceCpp
NULL

setClass (
  # Class name
  "Shallow3res",
  
  # Defining slot type
  representation (
    mut_mat = "matrix",
    wt = "numeric"
  ),
  
  # Initializing slots
  prototype = list(
    name = as.character(NULL),
    birth = as.Date(as.character(NULL))
  )
)

mut_count3 <- function(datapath) {
  # Make sure input is a data.table
  #load("./analysis/input2chromloop.Rda")
  datapath <- curated_lol 
  nb_uniq = length(unique(as.character(datapath$icgc_sample_id)))
  # This is the row names of init_mut_mat
  uniq_sample <- unique(as.character(datapath$icgc_sample_id))
  
  init_mut_mat <- matrix(0L, nrow = nb_uniq, ncol = 96) 
  init_wt <- rep(0, 96)
  
  # Call readfast
  # Check that no alleles are higher than 3 for mut and
  # Change all alternate alleles to their corresponding numbers
  # no sex chromosome
  # Call C++ function
  
  shallow3res <- new("Shallow3res", mut_mat = init_mut_mat, wt = init_wt)
  shallow3res = shallow_loop3(shallow3res)
}

#' @importFrom data.table fread data.table as.data.table
readfast <- function(datapath) {
  x <- data.table::fread(datapath, header = FALSE)
  invisible(x)
  
  # Handle factors if we are dealing with data in env.
}

#' @export
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

  datapath <- datapath[from_allele < 4 & mutated_to_allele < 4, .(icgc_sample_id, chromosome, 
                                            chromosome_start, 
                                            chromosome_end, mutated_from_allele,
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

#lol <- convSig:::readfast("~/Documents/GitHub/convSig-shallow/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
# hey <- fasta_process(lol$V1)
# hey <- fasta_process3(hello)
# load("~/Documents/GitHub/convSig-shallow/hello.Rda")
#class(lol$V1)
