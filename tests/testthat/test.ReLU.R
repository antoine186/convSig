context("Shallow ReLU Tests")

conv3 <- conv_create(32, 5, 3)
conv5 <- conv_create(512, 5, 5)

test_that("conv_create returns a well-formed result",
          {
            expect_equal(dim(conv3)[1], 32)
            expect_equal(dim(conv5)[1], 512)
            expect_equal(dim(conv3)[2], 5)
            expect_equal(dim(conv3)[2], 5)
          }
)

frag3 <- fragbase_indexer(3, 32)
frag5 <- fragbase_indexer(5, 512)

test_that("fragbase_indexer returns a well-formed result",
          {
            expect_equal(length(frag3), 3)
            expect_equal(length(frag3[[1]]), 4)
            expect_equal(length(frag3[[2]]), 2)
            expect_equal(length(frag3[[3]]), 4)
            expect_equal(length(frag5), 5)
            expect_equal(length(frag5[[1]]), 4)
            expect_equal(length(frag5[[2]]), 4)
            expect_equal(length(frag5[[3]]), 2)
            expect_equal(length(frag5[[4]]), 4)
            expect_equal(length(frag5[[5]]), 4)
          }
)

ten3 <- tencode(3, 32)
ten5 <- tencode(5, 512)

test_that("tencode returns a well-formed result",
          {
            expect_equal(dim(ten3)[1], 32)      
            expect_equal(dim(ten3)[2], 10)
            expect_equal(dim(ten5)[1], 512)      
            expect_equal(dim(ten5)[2], 18)
          }
)




