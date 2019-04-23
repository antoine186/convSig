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

test_csv <- 

test_that("ICGC2Mut works correctly with matrices and data.frames when data.loaded is set to TRUE",
          {
            
          }
)