context("Check ABC functions")
source("ABC.R")
source('test_helper.R')

qc = 0.6
qh = 0.4
n = 5
W <- prob_matrix(qc, qh, n)
obs_matrix <- matrix(c(66,87,25,22,4,
                     13,14,15,9,4,
                     0,4,4,9,1,
                     0,0,4,3,1,
                     0,0,0,1,1,
                     0,0,0,0,0), nrow = 6, byrow = TRUE)

test_that("All columns in probability matrix are equal to 1", {
  expect_equal(apply(W, 2, sum), rep(1, n+1))
})


test_that("Lower half of probability matrix is equal to 0", {
  expect_equal(W[lower.tri(W)], rep(0, n*(ceiling((n)/2))))
})

test_that("Probability distribution matches sampling distribution", {
  results <- abc_sim(obs_matrix, prior, 10, 100000)
  sum <- test_acc_params(table1, results1[[1]], results1[[2]])
  expect_true(sum < 1)
})
