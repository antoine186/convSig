context("ICGC Processing Tests")

test_that("base to number conversion works",
          {
            expect_equal(Base2Num("a"), 0)
            expect_equal(Base2Num("A"), 0)
            expect_equal(Base2Num("g"), 1)
            expect_equal(Base2Num("G"), 1)
            expect_equal(Base2Num("c"), 2)
            expect_equal(Base2Num("C"), 2)
            expect_equal(Base2Num("t"), 3)
            expect_equal(Base2Num("T"), 3)
            expect_equal(Base2Num("-"), 4)
            expect_equal(Base2Num("N"), 4)
            expect_equal(Base2Num(" "), 4)
            expect_equal(Base2Num("NA"), 4)
            expect_equal(Base2Num("AA"), 4)
            expect_equal(Base2Num("GGGA"), 4)
          }
)

valid_test_tsv <- data.table::fread("result_mutation_dataset.tsv")

test_that("icgc2mut is consistent for csv and tsv files",
          {
            expect_equal(icgc2mut("example_mutation_dataset.tsv",
                                  "GRCh37", "WGS"), valid_test_tsv)
            expect_equal(icgc2mut("example_mutation_dataset.csv",
                                  "GRCh37", "WGS"), valid_test_tsv)
          }
)

test_that("icgc2mut throws correct errors in erroneous argument scenarios",
          {
            expect_error(icgc2mut("example_mutation_dataset.tsv", "GRCh399999", "WGS"))
            expect_error(icgc2mut("example_mutation_dataset.tsv", "GRCh37", "WasdasdS"))
            expect_error(icgc2mut("example_mutati.tsv", "GRCh37", "WGS"))
          }
)

test_that("icgc2mut throws warning when non-mission critical arguments are missing", 
          {
            expect_warning(icgc2mut("example_mutation_dataset.tsv"))
            expect_warning(icgc2mut("example_mutation_dataset.tsv", 
                                    assembly = "GRCh37"))
            expect_warning(icgc2mut("example_mutation_dataset.tsv", Seq = "WGS"))
          }
)

example_mutation_dataset <- read.csv("example_mutation_dataset.csv")
example_mutation_dataset$assembly_version <- 
  as.character(example_mutation_dataset$assembly_version)

ins_ind <- c(1:20)
example_mutation_dataset$assembly_version[ins_ind] <- "world"
ins_ind <- c(21: 30)
example_mutation_dataset$assembly_version[ins_ind] <- "sugar"
ins_ind <- sample(31: 100)
example_mutation_dataset$assembly_version[ins_ind] <- "hello"

res_miss_assembly <- icgc2mut(example_mutation_dataset, assembly ="hello", Seq = "WGS")

test_that("icgc2mut behaves as expected when assembly argument is missing", 
          {
            expect_equal(icgc2mut(example_mutation_dataset, Seq = "WGS"), res_miss_assembly)
            expect_equal(icgc2mut(example_mutation_dataset), res_miss_assembly)
          }
)

example_mutation_dataset <- read.csv("example_mutation_dataset.csv")
example_mutation_dataset$sequencing_strategy <- 
  as.character(example_mutation_dataset$sequencing_strategy)

ins_ind <- c(1:20)
example_mutation_dataset$sequencing_strategy[ins_ind] <- "seq strat. 1"
ins_ind <- c(21: 30)
example_mutation_dataset$sequencing_strategy[ins_ind] <- "seq strat. 2"
ins_ind <- sample(31: 100)
example_mutation_dataset$sequencing_strategy[ins_ind] <- "seq strat. 3"

res_miss_seq <- icgc2mut(example_mutation_dataset, assembly ="GRCh37", Seq = "seq strat. 3")

test_that("icgc2mut behaves as expected when sequencing strat argument is missing", 
          {
            expect_equal(icgc2mut(example_mutation_dataset, assembly ="GRCh37"), res_miss_seq)
            expect_equal(icgc2mut(example_mutation_dataset), res_miss_seq)
          }
)

test_mat_frame <- read.csv("example_mutation_dataset.csv")
test_mat_frame$X <- NULL
test_mat_frame_mat <- as.matrix(test_mat_frame)
test_mat_frame_df <- as.data.frame(test_mat_frame)

test_that("icgc2mut works properly with matrices and data.frames 
          when data.loaded is set to TRUE",
          {
            expect_equal(icgc2mut(test_mat_frame_mat, "GRCh37", "WGS"), valid_test_tsv)
            expect_equal(icgc2mut(test_mat_frame_df, "GRCh37", "WGS"), valid_test_tsv)
          }
)

test_miss_from <- test_mat_frame
test_miss_to <- test_mat_frame
rm(test_mat_frame)
test_miss_from$mutated_from_allele <- NULL
test_miss_to$mutated_to_allele <- NULL

test_that("icgc2mut throws correct errors when mission critical arguments are missing", 
          {
            expect_error(icgc2mut(test_miss_from, "GRCh37", "WGS", data.loaded = TRUE))
            expect_error(icgc2mut(test_miss_to, "GRCh37", "WGS", data.loaded = TRUE))
          }
)

presort_input_same <- data.table::data.table(chromosome = c(1,3,2,1,3,3,2,1,2), 
                                             chromosome_start = c(11,24,9,12,14,16,21,32,4), 
                                             chromosome_end = c(11,24,9,12,14,16,21,32,4))

presort_input_notsame <- data.table::data.table(chromosome = c(1,3,2,1,3,3,2,1,2), 
                                                chromosome_start = c(11,24,9,12,14,16,21,32,4), 
                                                chromosome_end = c(12,11,17,10,16,8,32,31,9))

res_same <- data.table::data.table(chromosome = c(1,1,1,2,2,2,3,3,3), 
                                   chromosome_start = c(11,12,32,4,9,21,14,16,24),
                                   chromosome_end = c(11,12,32,4,9,21,14,16,24))

res_notsame <- data.table::data.table(chromosome = c(1,1,1,2,2,2,3,3,3),
                                      chromosome_start = c(11,12,32,4,9,21,14,16,24),
                                      chromosome_end = c(12,10,31,9,17,32,16,8,11))

test_that("sorting of rows works correctly",
          {
            expect_equal(icgc_sort(presort_input_same), res_same)
            expect_equal(icgc_sort(presort_input_notsame), res_notsame)
          }
)

input_same <- data.table::data.table(chromosome = c(1,3,2,1,3,3), 
                                     chromosome_start = c(11,24,9,12,14,16), 
                                     chromosome_end = c(11,24,9,12,14,16))

input_notsame <- data.table::data.table(chromosome = c(1,3,2,1,3,3), 
                                                chromosome_start = c(11,24,9,12,14,16), 
                                                chromosome_end = c(11,24,17,12,15,16))

res_same <- data.table(chromosome = c(1,3,2,1,3,3), 
                                   chromosome_start = c(11,24,9,12,14,16),
                                   chromosome_end = c(11,24,9,12,14,16))

res_notsame <- data.table::data.table(chromosome = c(1,3,1,3),
                                      chromosome_start = c(11,24,12,16),
                                      chromosome_end = c(11,24,12,16))

input_allnotsame <- data.table::data.table(chromosome = c(1,3,2), 
                                     chromosome_start = c(11,24,9), 
                                     chromosome_end = c(12,25,10))

test_that("that the removal of non single nucleotide changes works properly",
          {
            expect_equal(icgc_snp(input_same), res_same)
            expect_equal(icgc_snp(input_notsame), res_notsame)
            expect_error(icgc_snp(input_allnotsame))
          }
)

dup_input <- data.table::data.table(x=c(1,1,1,1), y=c(1,2,3,1))
dup_res <- data.table::data.table(x=c(1,1,1), y=c(1,2,3))

test_that("removal of duplicated rows works correctly",
          {
            expect_equal(icgc_dedup(dup_input), dup_res)
          }
)

# Test throwing errors when input is malformed
# When are not data.frame or matrix
# When missing critical columns

res <- icgc2mut("example_mutation_dataset.tsv", "GRCh37", "WGS")
test_res_mat <- as.matrix(res)
test_res_df <- as.data.frame(res)

res_rmnonSNP <- icgc_curate(res, remove.nonSNP = TRUE)

test_that("icgc_curate returns consistent results for both data.frames and matrices",
          {
            expect_equal(icgc_curate(test_res_mat, remove.nonSNP = TRUE), res_rmnonSNP)
            expect_equal(icgc_curate(test_res_df, remove.nonSNP = TRUE), res_rmnonSNP)
          }  
)

res_nochrom <- res
res_nostart <- res
res_noend <- res
res_nochrom$chromosome <- NULL
res_nostart$chromosome <- NULL
res_noend$chromosome <- NULL

test_that("icgc_curate throws appropriate errors when the input is malformed",
          {
            expect_error(icgc_curate("incorrect input", remove.nonSNP = TRUE))
            expect_error(icgc_curate(res_nochrom, remove.nonSNP = TRUE))
            expect_error(icgc_curate(res_nostart, remove.nonSNP = TRUE))
            expect_error(icgc_curate(res_noend, remove.nonSNP = TRUE))
          }
)


