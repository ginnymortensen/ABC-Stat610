# -----------------------------------
#
#title: STAT610 FA2023 Project
#author: Ben Iovino
#date: 12/3/23
#
# -----------------------------------


prior <- function() {
  return(runif(1))
}


prob_matrix <- function(qc, qh, n) {
  #' Returns a nxn matrix of probabilities given two values
  #' 
  #' @param qc probability 1
  #' @param qh probability 2
  #' @param n number of rows and columns
  #' @return matrix of probabilities
  
  # Initialize matrix with row/col starting at 0
  W = matrix(0, nrow = n+1, ncol = n+1)
  
  # First row is qc^s, s = 0, 1, ..., n
  W[1,] = qc^(0:n)
  
  # Subsequent rows
  for (j in 2:(n+1)) {
    for (i in j:(n+1)) {
      if (i == j) {
        W[j,j] <- 1 - sum(W[1:(j-1), j])
      }
      else {
        W[j, i] <- choose(i-1,j-1) * W[j,j] * (qc*qh^(j-1))^(i-j)
      }
    }
  }
  
  return(W)
  
}


sample_from_matrix_distribution <- function(prob_matrix, obs_matrix) {
  #' Returns a matrix of sampled data given a probability matrix and observed data
  #' 
  #' @param prob_matrix probability matrix (nxn)
  #' @param obs_matrix observed data (nxn)
  #' @return matrix of sampled data

  prob_matrix <- prob_matrix[, -1]  # 1st column not relevant
  sampled_data <- matrix(0, nrow = nrow(obs_matrix), ncol = ncol(obs_matrix))
  
  for (c in 1:ncol(sampled_data)) {
    sim_vector <- rmultinom(1, sum(obs_matrix[,c]), prob = prob_matrix[,c])
    sampled_data[,c] <- sim_vector
  }
  
  return(sampled_data)
  
}


abc_sim <- function(obs_data, prior, tolerance, n_steps) {
  #' Returns a list of two vectors of accepted parameters from data simulation
  #' 
  #' @param obs_data observed data (nxn)
  #' @param prior function that samples qc and qh from prior distribution
  #' @param tolerance acceptance tolerance for simulated data
  #' @param n_steps number of steps to run ABC
  #' @return list of two vectors

  qc_accepted <- c()
  qh_accepted <- c()
  for (i in 1:n_steps) {  # Simulate data n_step times
    qc <- prior()
    qh <- prior()
    n = ncol(obs_data)
    W <- prob_matrix(qc, qh, n)
    sim_data <- sample_from_matrix_distribution(W, obs_data)
    if (norm(obs_data - sim_data, "F") <= tolerance) {  # F-norm against tolerance
      qc_accepted <- append(qc_accepted, qc)
      qh_accepted <- append(qh_accepted, qh)
    }
  }
  
  return(list(qh_accepted, qc_accepted))
}


plot_results <- function(results1, results2, title, legend1, legend2) {
  #' Plots results from ABC simulation
  #' 
  #' @param results1 list of two vectors of accepted parameters from data simulation
  #' @param results2 list of two vectors of accepted parameters from data simulation
  #' @param title title of plot
  #' @param legend1 legend for first vector
  #' @param legend2 legend for second vector
  
  plot(results1[[1]], results1[[2]], col='red', xlim=c(0,1), ylim=c(0,1),
       xlab="qH", ylab="qC", main=title)
  points(results2[[1]], results2[[2]], col='blue')
  legend("topleft", legend=c(legend1, legend2), col=c("red", "blue"), pch=1)
}


# Read csv files
table1 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table1.csv"))
table2 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table2.csv"))
table3 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table3.csv"))
table4 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table4.csv"))

# Run ABC sim for table s2 and plot results
tolerance = 20
n_steps = 100000
results1 <- abc_sim(table1, prior, tolerance, n_steps)
results2 <- abc_sim(table2, prior, tolerance, n_steps)
plot_results(results1, results2,
             "ABC Posterior Distributions for Table S2", "1977-78", "1980-81")

# Run ABC sim for table s3 and plot results
tolerance = 10
n_steps = 100000
results1 <- abc_sim(table3, prior, tolerance, n_steps)
results2 <- abc_sim(table4, prior, tolerance, n_steps)
plot_results(results1, results2,
             "ABC Posterior Distributions for Table S3", "1975-76", "1978-79")

W <- prob_matrix(0.6, 0.4, 5)
print(W)

