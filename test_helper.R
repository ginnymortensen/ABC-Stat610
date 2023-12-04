# -----------------------------------
#
#title: STAT610 FA2023 Project
#author: Ben Iovino
#date: 12/4/23
#
# -----------------------------------


sample_table_sum <- function(samp, W) {
  #' Returns the sum of the absolute difference between each index in the
  #'sampled table and the probability matrix.
  #' 
  #' @param samp sampled table
  #' @param W probability matrix
  #' @return sum of absolute difference
  
  for (i in 1:ncol(samp)) {
    samp[,i] = samp[,i]/sum(samp[,i])
    samp[,i] = abs(samp[,i] - W[,i])
  }

  return(sum(samp))
}


test_acc_params <- function(obs, qc, qh) {
  #' Generates probability matrices for each value in qc and qh. Compares the
  #' sampled data to the observed data and returns the average sum between
  #' the absolute differences of the two.
  #' 
  #' @param obs observed data
  #' @param qc vector of accepted parameter values
  #' @param qh vector of accepted parameter values
  #' @return average sum of absolute differences between sampled and observed data
  
  n = ncol(obs)
  sum = 0
  for (i in 1:length(qc)) {
    W = prob_matrix(qc[i], qh[i], n)
    sim_data = sample_from_matrix_distribution(W, obs)
    sum = sum + sample_table_sum(sim_data, W[,-1])
  }
  
  return(sum/length(qc))
  
}
