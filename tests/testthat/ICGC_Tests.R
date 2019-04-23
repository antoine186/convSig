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

# Should remove those later
input_test_tsv <- data.table::fread("./inst/extdata/example_mutation_dataset.tsv")
valid_test_tsv <- data.table::fread("./inst/extdata/result_mutation_dataset.tsv")

test_that("ICGC2Mut is consistent for csv and tsv files",
          {
            expect_equal(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", "GRCh37", "WGS"), valid_test_tsv)
            expect_equal(ICGC2Mut("./inst/extdata/example_mutation_dataset.csv", "GRCh37", "WGS"), valid_test_tsv)
          }
)

test_that("ICGC2Mut throws correct errors in erroneous argument scenarios",
          {
            expect_equal(ICGC2Mut("./inst/extdata/example_mutation_dataset.tsv", "GRCh37", "WGS"), valid_test_tsv)
          }
)

