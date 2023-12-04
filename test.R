context("Check ABC functions")
source("ABC.R")

qc = 0.6
qh = 0.4
n = 5
W <- prob_matrix(qc, qh, n)

test_that("All columns in probability matrix are equal to 1", {
  expect_equal(apply(W, 2, sum), rep(1, n+1))
})


test_that("Lower half of probability matrix is equal to 0", {
  expect_equal(W[lower.tri(W)], rep(0, n*(ceiling((n)/2))))
})

