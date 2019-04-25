context("ICGC2Mut Tests")

library(convSig)

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

valid_test_tsv <- data.table::fread("./inst/extdata/result_mutation_dataset.tsv")

test_that("ICGC2Mut is consistent for csv and tsv files",
          {
            expect_equal(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", "GRCh37", "WGS"), valid_test_tsv)
            expect_equal(ICGC2Mut("./inst/extdata/example_mutation_dataset.csv", "GRCh37", "WGS"), valid_test_tsv)
          }
)

test_that("ICGC2Mut throws correct errors in erroneous argument scenarios",
          {
            expect_error(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", "GRCh399999", "WGS"))
            expect_error(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", "GRCh37", "WasdasdS"))
            expect_error(ICGC2Mut("./inst/extdata/example_mutati.tsv", "GRCh37", "WGS"))
          }
)

test_that("ICGC2Mut throws warning when non-mission critical arguments are missing", 
          {
            expect_warning(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv"))
            expect_warning(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", assembly = "GRCh37"))
            expect_warning(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", Seq = "WGS"))
          }
)

test_mat_frame <- read.csv("./inst/extdata/example_mutation_dataset.csv")
test_mat_frame$X <- NULL
test_mat_frame_mat <- as.matrix(test_mat_frame)
test_mat_frame_df <- as.data.frame(test_mat_frame)

test_that("ICGC2Mut works properly with matrices and data.frames when data.loaded is set to TRUE",
          {
            expect_equal(ICGC2Mut(test_mat_frame_mat, "GRCh37", "WGS", data.loaded = TRUE), valid_test_tsv)
            expect_equal(ICGC2Mut(test_mat_frame_df, "GRCh37", "WGS", data.loaded = TRUE), valid_test_tsv)
          }
)

rm(valid_test_tsv)
rm(test_mat_frame_mat)
rm(test_mat_frame_df)

test_miss_from <- test_mat_frame
test_miss_to <- test_mat_frame
rm(test_mat_frame)
test_miss_from$mutated_from_allele <- NULL
test_miss_to$mutated_to_allele <- NULL

test_that("ICGC2Mut throws correct errors when mission critical arguments are missing", 
          {
            expect_error(ICGC2Mut(test_miss_from, "GRCh37", "WGS", data.loaded = TRUE))
            expect_error(ICGC2Mut(test_miss_to, "GRCh37", "WGS", data.loaded = TRUE))
          }
)

rm(test_miss_from)
rm(test_miss_to)

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
            expect_equal(ICGC_sort(presort_input_same), res_same)
            expect_equal(ICGC_sort(presort_input_notsame), res_notsame)
          }
)

rm(presort_input_same)
rm(presort_input_notsame)
rm(res_same)
rm(res_notsame)

dup_input <- data.table::data.table(x=c(1,1,1,1), y=c(1,2,3,1))
dup_res <- data.table::data.table(x=c(1,1,1), y=c(1,2,3))

test_that("removal of duplicated rows works correctly",
          {
            expect_equal(ICGC_dedup(dup_input), dup_res)
          }
)

rm(dup_input)
rm(dup_res)


