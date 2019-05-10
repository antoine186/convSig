context("Shallow Trinucleotide Window Tests")

interdata <- icgc2mut("example_mutation_dataset.tsv",
                      "GRCh37", "WGS")
interdata <- icgc_curate(interdata)

rand_insert <- c(34, 18, 23, 6)
interdata$chromosome[rand_insert] = "X"
rand_insert <- c(27, 7, 21, 36)
interdata$chromosome[rand_insert] = "Y"

rand_insert <- c(17)
interdata$mutated_from_allele[rand_insert] = "N"
rand_insert <- c(5)
interdata$mutated_to_allele[rand_insert] = "N"

from_allele <- lapply(interdata[, mutated_from_allele], Base2Num)
to_allele <- lapply(interdata[, mutated_to_allele], Base2Num)

rm_from <- which(from_allele == 4)
rm_to <- which(to_allele == 4)

res_nosexorn <- interdata[-rm_from,]
res_nosexorn <- res_nosexorn[-rm_to,]

res_nosexorn <- res_nosexorn[!chromosome == "X"]
res_nosexorn <- res_nosexorn[!chromosome == "Y"]

res_nosexorn$mutated_to_allele <- NULL

test_that("skipping sex chromosomes and padding \"N\" bases works",
          {
            x <- mut_process3(interdata)
            x$mutated_to_allele <- NULL
            expect_equal(x, res_nosexorn)
          }
)

testdata <- icgc2mut("example_mutation_dataset.tsv",
                      "GRCh37", "WGS")
testdata <- icgc_curate(testdata)

interdata <- testdata

from_allele <- lapply(interdata[, mutated_from_allele], Base2Num)
to_allele <- lapply(interdata[, mutated_to_allele], Base2Num)

from_allele <- as.numeric(from_allele)
to_allele <- as.numeric(to_allele)

for (i in 1:dim(interdata)[1]) {
  if (from_allele[i] < 2) {
    to_allele[i] = 3 - to_allele[i]
  }
  if (to_allele[i] > 2) {
    to_allele[i] = 2
  }
}

interdata$mutated_to_allele <- to_allele

test_that("alternate allele transformation works",
          {
            x <- mut_process3(testdata)
            expect_equal(x, interdata)
          }
)
