# Clean R Script
rm(list = ls())

# Load functions
source("functions.R")

# Load Required Libraries
library(sandwich)
library(compiler)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(plm)
library(tictoc)

# Set Seed for Reproducibility
set.seed(1)

# Constants
Rep <- 10000
N <- 50
T <- 50
config_mat <- matrix(c(
  1, 0, 0, 1, 0.7, 0,
  1, 0.7, 0.7, 2, 0, 0,
  2, 0.7, 0.7
), nrow = 3, ncol = 6)

# Start Timer
tic()

# Determine the number of cores available and use all but 4
total_cores <- detectCores()
cores_to_use <- total_cores - 4
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# Initialize Result Matrix
save_result <- matrix(0, nrow = 6, ncol = 7)

# Main Simulation Loop
for (j in seq_len(ncol(config_mat))) {
  DGP <- config_mat[1, j]
  rho <- config_mat[2, j]
  kappa <- config_mat[3, j]

  res <- foreach(k = 1:Rep,
                 .combine = 'rbind',
                 .errorhandling = 'remove',
                 .packages = c("sandwich", "plm")) %dorng% {
                   # Generate Random Variables
                   gen <- generatei(N, T, DGP, rho, kappa)
                   data <- gen$data

                   # Formula for regression
                   formula <- vY ~ vX - 1

                   # Baseline estimate
                   est <- estimator_dci(formula, data)
                   std_errors <- sqrt(est$beta_var)
                   Save <- est$beta
                   Cov <- (2 <= Save + 1.96 * std_errors / sqrt(N)) &
                     (2 >= Save - 1.96 * std_errors / sqrt(N))
                   Wid <- 2 * 1.96 * std_errors / sqrt(N)

                   cbind(Save, Cov, Wid)
                 }

  # Store Results
  save_result[j, ] <- c(
    DGP,
    rho,
    kappa,
    mean(res[, 1]) - 2,                      # Bias
    mean((res[, 1] - mean(res[, 1]))^2),     # Variance
    mean(res[, 2]),                          # Coverage
    mean(res[, 3])                           # Width
  )
}

# Convert results to data frame
results <- data.frame(save_result)

# Set column names
colnames(results) <- c("DGP", "rho", "kappa", "Bias", "Var", "Cov", "Wid")

# Save results to file
filename <- paste("results_simi", "_N", N, "_T", T, ".RData", sep = "")
save(results, file = filename)

# Display results
results

# End Timer
toc()
