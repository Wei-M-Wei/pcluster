rm(list = ls())

# Load Required Libraries
library(sandwich)
library(gtools)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(tictoc)
library(pcluster)

# Set Seed for Reproducibility
set.seed(1)

# Generate data

N <- 20
T <- 20
E <- matrix( rnorm(N*T), N, T)
U <- matrix( rnorm(N*T), N, T)
E[, 1] <- rnorm(N)
U[, 1] <- rnorm(N)


A <- rgamma(N, 1)
B <- rgamma(N, 1)

A_rep <- matrix(rep(A, T), N, T)
B_rep <- t(matrix(rep(B, N), T, N))


F <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 10)
H <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 5)


X <- H + U
Y <- X + F + E

id_code <- numeric(N * T)
for (i in 1:N) {
  for (t in 1:T) {
    id_code[((t - 1) * N + i):(t * N)] <- 10+i
  }
}

time_code <- numeric(N * T)
for (i in 1:N) {
  for (t in 1:T) {
    time_code[((t - 1) * N + 1):(t * N)] <- 10 + t
  }
}

vY <- c(Y)
vX <- c(X)

data <- data.frame(id_code, time_code, vY, vX)


formula <- vY ~ vX - 1
index = c("id_code", "time_code")
init <- 300

# Baseline estimate
help("estimator_dc")

est <- estimator_dc(formula, data, index, init = init)
ols <- est[["res"]]
G <- est[["G"]]
C <- est[["C"]]
summary_correct = est$summary_table
print(summary_correct)
coef_estimate = summary_correct$coefficients$Estimate

# We use the heteroskedasticity autocorrelation consistent standard errors clustered at the level of each unit together with correction for the degrees of freedom, see our paper page 20.
std_error_corrected = summary_correct$coefficients$`Std. Error corrected`


# Cross-fitted estimate
est_CF <- estimator_dc(formula, data, index, CF = TRUE, init = init)
ols_CF <- est_CF[["res"]]
coeff_CF <- ols_CF$coefficients
summary_correct_CF = est_CF$summary_table
print(summary_correct_CF)
coef_estimate_CF = summary_correct_CF$coefficients$Estimate
std_error_corrected_CF = summary_correct_CF$coefficients$`Std. Error corrected`
