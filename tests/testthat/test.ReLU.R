context("Shallow ReLU Tests")

K = 5
poss3 = (3 * 4) - 2
feat3 <- matrix(runif(poss3*K), nrow = poss3, ncol = K)
poss5 = (5 * 4) - 2
feat5 <- matrix(runif(poss5*K), nrow = poss5, ncol = K)
complexfeat3 <- array(runif(3 * 4 * K), dim = c(3, 4, K))
complexfeat5 <- array(runif(5 * 4 * K), dim = c(5, 4, K))
conv3 <- conv_create(32, 5, 3, feat3)
conv5 <- conv_create(512, 5, 5, feat5)
complexconv3 <- complexconv_create(32, 5, 3, complexfeat3, 2)
complexconv5 <- complexconv_create(512, 5, 5, complexfeat5, 3)

test_that("conv_create returns a well-formed result",
          {
            expect_equal(dim(conv3)[1], 32)
            expect_equal(dim(conv5)[1], 512)
            expect_equal(dim(conv3)[2], 5)
            expect_equal(dim(conv3)[2], 5)
          }
)

test_that("complexconv_create returns a well-formed result",
          {
            expect_equal(dim(complexconv3)[1], 32)
            expect_equal(dim(complexconv5)[1], 512)
            expect_equal(dim(complexconv3)[2], 5)
            expect_equal(dim(complexconv3)[2], 5)
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

ar1 <- array(runif(10*5), dim = c(10,1,5))
ar2 <- array(runif(10*5), dim = c(10,5,1))
ar3 <- array(runif(10*5), dim = c(1,10,5))

test_that("three_recycle produces the correct format at any designated axis",
          {
            expect_equal(dim(three_recycle(ar1, 2, 20)), c(10, 20, 5))
            expect_equal(dim(three_recycle(ar2, 3, 20)), c(10, 5, 20))   
            expect_equal(dim(three_recycle(ar3, 1, 20)), c(20, 10, 5))   
          }
)

ar1 <- array(runif(10*5*3), dim = c(10,3,5,1))
ar2 <- array(runif(10*5*3), dim = c(10,5,1,3))
ar3 <- array(runif(10*5*3), dim = c(3,1,10,5))
ar4 <- array(runif(10*5*3), dim = c(1,3,10,5))

test_that("four_recycle produces the correct format at any designated axis",
          {
            expect_equal(dim(four_recycle(ar1, 4, 20)), c(10,3,5,20))
            expect_equal(dim(four_recycle(ar2, 3, 20)), c(10,5,20,3))
            expect_equal(dim(four_recycle(ar3, 2, 20)), c(3,20,10,5))
            expect_equal(dim(four_recycle(ar4, 1, 20)), c(20,3,10,5))
          }
)

ar <- array(runif(10*3*5), dim = c(10,3,5))
ar_summed1 <- apply(ar,c(1,3),sum)
ar_summed2 <- apply(ar,c(1,2),sum)

test_that("three_colsum produces the correct result",
          {
            expect_equal(array(three_colsum(ar, 2), dim = c(10,5)), ar_summed1)
            expect_equal(array(three_colsum(ar, 3), dim = c(10,3)), ar_summed2)
          }
)

ar <- array(runif(10*3*5*8), dim = c(10,3,5,8))
ar_summed1 <- apply(ar,c(1,2,3),sum)
ar_summed2 <- apply(ar,c(1,2,4),sum)
ar_summed3 <- apply(ar,c(1,3,4),sum)

test_that("four_colsum produces the correct format",
          {
            expect_equal(array(four_colsum(ar, 4), dim = c(10,3,5)), ar_summed1)
            expect_equal(array(four_colsum(ar, 3), dim = c(10,3,8)), ar_summed2)
            expect_equal(array(four_colsum(ar, 2), dim = c(10,5,8)), ar_summed3)
          }
)

ar <- array(runif(10*3*5*8), dim = c(10,3,5,8,12))
ar_summed1 <- apply(ar,c(1,2,4,5),sum)

test_that("five_colsum produces the correct format",
          {
            expect_equal(array(five_colsum(ar, 3), dim = c(10,3,8,12)), ar_summed1)
          }
)

X <- array(runif(10*32*3), dim = c(10,32,3))
Z <- array(runif(10*32*3*5), dim = c(10,32,3,5))
type <- fragbase_indexer(5, 32)
mid = 2
N = 32

inter_X = list_indexer(X, type, mid, N, K)
inter_Z = list_indexer(Z, type, mid, N, K, X = FALSE)

test_that("list_indexer works properly for 3 bases for both input types",
          {
            expect_equal(dim(inter_X), c(10,2,16,3))
            expect_equal(dim(inter_Z), c(10,2,16,3, 5))
          }
)

X <- array(runif(10*512*3), dim = c(10,512,3))
Z <- array(runif(10*512*3*5), dim = c(10,512,3,5))
type <- fragbase_indexer(5, 512)
mid = 3
N = 512

inter_X = list_indexer(X, type, mid, N, K)
inter_Z = list_indexer(Z, type, mid, N, K, X = FALSE)

test_that("list_indexer works properly for 5 bases for both input types",
          {
            expect_equal(dim(inter_X), c(10,2,256,3))
            expect_equal(dim(inter_Z), c(10,2,256,3, 5))
          }
)



