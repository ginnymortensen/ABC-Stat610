# -----------------------------------
#
#title: STAT610 FA2023 Project
#author: Ben Iovino
#date: 12/?/23
#
# -----------------------------------

# Read csv files
table1 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table1.csv"))
table2 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table2.csv"))
table3 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table3.csv"))
table4 = as.matrix(read.csv("/home/ben/Code/STAT610/Project/data/table4.csv"))

table1

# Randomly sample qc and qh from uniform distribution
qc = 0.6
qh = 0.4
n = 5

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
        W[j, i] <- choose(i-1,j-1) * W[j,j] * (qc*qh^j-1)^((i-1)-(j-1))
      }
    }
  }
  
  return(W)
  
}

W <- prob_matrix(qc, qh, n)
print(W)

F_norm <- norm(W, type="F")




