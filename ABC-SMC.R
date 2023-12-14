source("ABC.R")


kernel <- function(population, particle) {
  #' Returns the perturbed particle 
  #' @param population previous population
  #' @param particle particle to be perturbed
  max_val <- max(population)
  min_val <- min(population)
  sigma <- 0.5*(max_val - min_val)
  return(particle + runif(1, -sigma, sigma))
}


# Checks whether the particle lies within the support of prior 
prior_pdf <- function(x) {
  return(x >= 0 && x <= 1)
}


get_param_sample <- function(population, prior_pdf) {
  #' Returns the perturbed particle lying within the support of prior
  #' 
  #' @param population previous population
  #' @param prior_pdf function to verify that the particle lies within the support of prior
  while (TRUE) {
    particle <- sample(population, 1)
    perturbed <- kernel(population, particle)
    if (prior_pdf(perturbed)) {
      return(perturbed)
    }
  }
}


abc_smc_final <- function(obs_data1, prior, tolerance, n_steps) {
  #' Returns a list of two vectors of accepted parameters from data simulation
  #' 
  #' @param obs_data observed data (nxn)
  #' @param prior function that samples qc and qh from prior distribution
  #' @param tolerance list of acceptance tolerance for simulated data
  #' @param n_steps number of particles to accept for each tolerance value
  #' @return list of two vectors
  n1 = ncol(obs_data1)
  for (i in 1:(length(tolerance)-1)) {
    accepted_qc1 <- c()
    accepted_qh1 <- c()
    acc_count <- 0
    while (acc_count < n_steps) {
      if (i == 1) {
        qc1 <- prior()
        qh1 <- prior()
      }
      else {
        qc1 <- get_param_sample(qc1_population, prior_pdf)
        qh1 <- get_param_sample(qh1_population, prior_pdf)
      }
      W1 <- prob_matrix(qc1, qh1, n1)
      sim_data1 <- sample_from_matrix_distribution(W1, obs_data1)
      if (norm(obs_data1 - sim_data1, "F") <= tolerance[i]) {
        accepted_qc1 <- append(accepted_qc1, qc1)
        accepted_qh1 <- append(accepted_qh1, qc1)
        acc_count <- acc_count + 1
      }
    }
    qc1_population <- accepted_qc1
    qh1_population <- accepted_qh1
  }
  return(list(accepted_qh1, accepted_qc1))
}


# Read csv files
table1 = as.matrix(read.csv("D:\\sanjana\\table1.csv"))
table2 = as.matrix(read.csv("D:\\sanjana\\table2.csv"))
table3 = as.matrix(read.csv("D:\\sanjana\\table3.csv"))
table4 = as.matrix(read.csv("D:\\sanjana\\table4.csv"))


# Run ABC SMC for table s2 and plot results
tolerance = c(100, 80, 50, 30, 20, 15, 13, 12)
n_steps = 1000
results1 <- abc_smc_final(table1, prior, tolerance, n_steps)
results2 <- abc_smc_final(table2, prior, tolerance, n_steps)
plot_results(results1, results2,
             "ABC-SMC Posterior Distributions for Table S2", "1977-78", "1980-81")

# Run ABC SMC for table s3 and plot results
tolerance = c(40, 20, 15, 10, 8, 6, 5)
n_steps = 1000
results1 <- abc_smc_final(table3, prior, tolerance, n_steps)
results2 <- abc_smc_final(table4, prior, tolerance, n_steps)
plot_results(results1, results2,
             "ABC-SMC Posterior Distributions for Table S3", "1975-76", "1978-79")