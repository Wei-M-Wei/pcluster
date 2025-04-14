rm(list = ls())

# Load functions
source("functions.R")

# Load Required Libraries
library(sandwich)
library(gtools)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(plm)
library(pcluster)
library(tictoc)

# Set Seed for Reproducibility
set.seed(1)

# Constants
DGP <- 1
rho <- 0.7
kappa <- 0.0
Rep <- 10
Nlist <- c(50)
Tlist <- c(10, 20, 30, 40, 50)

# Start Timer
tic()

# Determine the number of cores available and use all but 4
total_cores <- detectCores()
cores_to_use <- total_cores - 4
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# Initialize Result Matrix
save_result <- matrix(0, nrow = length(Nlist) * length(Tlist), ncol = 48)

# Main Simulation Loop
for (n in seq_along(Nlist)) {
  for (t in seq_along(Tlist)) {
    N <- Nlist[n]
    T <- Tlist[t]

    res <- foreach(k = 1:Rep, .combine = 'rbind', .errorhandling = 'remove', .packages = c("gtools", "sandwich", "plm", "pcluster")) %dorng% {
      # Generate Random Variables
      data <- generate(N, T, DGP, rho, kappa)

      # Formula for regression
      formula <- vY ~ vX - 1

      # Baseline estimate
      est <- estimator_dc(formula, data, index = c("id", "time"))
      ols <- est[["res"]]
      G <- est[["G"]]
      C <- est[["C"]]
      std_errors <- est$summary_table$coefficients$`Std. Error corrected`

      Save <- ols$coefficients
      Cov <- (1 <= ols$coefficients + 1.96 * sqrt(N * T / ((N - G) * (T - C))) * std_errors) &
        (1 >= ols$coefficients - 1.96 * sqrt(N * T / ((N - G) * (T - C))) * std_errors)
      Wid <- 2 * 1.96 * sqrt(N * T / ((N - G) * (T - C))) * std_errors

      # Cross-fitted estimate
      est_CF <- estimator_dc(formula, data, CF = TRUE, index = c("id", "time"))

      ols_CF <- est_CF[["res"]]
      std_errors_CF <- est_CF$summary_table$coefficients$`Std. Error corrected`

      Save_CF <- ols_CF$coefficients
      Cov_CF <- (1 <= ols_CF$coefficients + 1.96 * std_errors_CF) &
        (1 >= ols_CF$coefficients - 1.96 * std_errors_CF)
      Wid_CF <- 2 * 1.96 * std_errors_CF

      # TWFE estimate
      ols_twfe <- plm(vY ~ vX, data = data, index = c("id", "time"), model = "within", effect = "twoways")
      std_errors_twfe <- sqrt(vcovHC(ols_twfe, type = "HC0", method = "arellano"))

      Save_twfe <- ols_twfe$coefficients
      Cov_twfe <- (1 <= ols_twfe$coefficients + 1.96 * std_errors_twfe) &
        (1 >= ols_twfe$coefficients - 1.96 * std_errors_twfe)
      Wid_twfe <- 2 * 1.96 * std_errors_twfe

      # Bai's estimate with R = T^{1/2}
      formula_vars <- all.vars(formula)
      dependent_var <- formula_vars[1]
      independent_vars <- formula_vars[-1]

      Y <- reshape_to_matrix(data, "id", "time", dependent_var)
      X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
      X <- do.call(cbind, X_list)
      R <- floor(T^(1/2))
      est_Bai <- LS.factor(X, Y, R = R)
      ols_Bai <- est_Bai[["beta"]]
      std_errors_Bai <- est_Bai[["hac.se"]]

      Save_Bai <- ols_Bai
      Cov_Bai <- (1 <= ols_Bai + 1.96 * std_errors_Bai) &
        (1 >= ols_Bai - 1.96 * std_errors_Bai)
      Wid_Bai <- 2 * 1.96 * std_errors_Bai

      # Estimate GFE
      R <- 5
      ls.output <- LS.factor(X, Y, R = R)
      X <- array(X, dim = c(nrow(X), ncol(X), 1))
      L <- ls.output$svd$u[, 1:5]
      F <- ls.output$svd$v[, 1:5]
      est_GFE <- hnn.est(Y, array(X, dim = c(N, T, 1)), L, F)
      ols_GFE <- est_GFE[["beta.hnn"]]
      std_errors_GFE <- est_GFE[["cluster.se"]]

      Save_GFE <- ols_GFE
      Cov_GFE <- (1 <= ols_GFE + 1.96 * std_errors_GFE) &
        (1 >= ols_GFE - 1.96 * std_errors_GFE)
      Wid_GFE <- 2 * 1.96 * std_errors_GFE

      # FA estimate
      est_FA <- FA(formula, data)
      ols_FA <- est_FA[["res"]]
      std_errors_FA <- sqrt(vcovHC(ols_FA, type = "HC0", method = "arellano"))

      Save_FA <- ols_FA$coefficients
      Cov_FA <- (1 <= ols_FA$coefficients + 1.96 * std_errors_FA) &
        (1 >= ols_FA$coefficients - 1.96 * std_errors_FA)
      Wid_FA <- 2 * 1.96 * std_errors_FA

      # One-way estimate
      est_BLM1 <- est_BLM_one_way(formula, data)
      ols_BLM1 <- est_BLM1[["res"]]
      std_errors_BLM1 <- sqrt(vcovHC(ols_BLM1, type = "HC0", method = "arellano")) * sqrt((N * T) / (N * T - est_BLM1[["df"]]))

      Save_BLM1 <- ols_BLM1$coefficients
      Cov_BLM1 <- (1 <= ols_BLM1$coefficients + 1.96 * std_errors_BLM1) &
        (1 >= ols_BLM1$coefficients - 1.96 * std_errors_BLM1)
      Wid_BLM1 <- 2 * 1.96 * std_errors_BLM1

      # Two-way estimate
      est_BLM2 <- est_BLM_two_way(formula, data)
      ols_BLM2 <- est_BLM2[["res"]]
      std_errors_BLM2 <- sqrt(vcovHC(ols_BLM2, type = "HC0", method = "arellano")) * sqrt((N * T) / (N * T - est_BLM2[["df"]]))

      Save_BLM2 <- ols_BLM2$coefficients
      Cov_BLM2 <- (1 <= ols_BLM2$coefficients + 1.96 * std_errors_BLM2) &
        (1 >= ols_BLM2$coefficients - 1.96 * std_errors_BLM2)
      Wid_BLM2 <- 2 * 1.96 * std_errors_BLM2

      # CCE estimator
      ccep <- pcce(vY ~ vX, data = data, model = "p")
      std_errors_cce <- sqrt(vcovHC(ccep, type = "HC0", method = "arellano"))

      Save_cce <- ccep$coefficients
      Cov_cce <- (1 <= ccep$coefficients + 1.96 * std_errors_cce) &
        (1 >= ccep$coefficients - 1.96 * std_errors_cce)
      Wid_cce <- 2 * 1.96 * std_errors_cce

      # Bai's estimate with R = T^{1/3}
      R <- floor(T^(1/3))
      est_Bai13 <- LS.factor(X, Y, R = R)
      ols_Bai13 <- est_Bai13[["beta"]]
      std_errors_Bai13 <- est_Bai13[["hac.se"]]

      Save_Bai13 <- ols_Bai13
      Cov_Bai13 <- (1 <= ols_Bai13 + 1.96 * std_errors_Bai13) &
        (1 >= ols_Bai13 - 1.96 * std_errors_Bai13)
      Wid_Bai13 <- 2 * 1.96 * std_errors_Bai13

      # Bai's estimate with R = T^{1/4}
      R <- floor(T^(1/4))
      est_Bai14 <- LS.factor(X, Y, R = R)
      ols_Bai14 <- est_Bai14[["beta"]]
      std_errors_Bai14 <- est_Bai14[["hac.se"]]

      Save_Bai14 <- ols_Bai14
      Cov_Bai14 <- (1 <= ols_Bai14 + 1.96 * std_errors_Bai14) &
        (1 >= ols_Bai14 - 1.96 * std_errors_Bai14)
      Wid_Bai14 <- 2 * 1.96 * std_errors_Bai14

      cbind(Save, Cov, Wid, G, C, Save_CF, Cov_CF, Wid_CF, Save_twfe, Cov_twfe, Wid_twfe, Save_Bai, Cov_Bai, Wid_Bai, Save_GFE, Cov_GFE, Wid_GFE, Save_FA, Cov_FA, Wid_FA, Save_BLM1, Cov_BLM1, Wid_BLM1, Save_BLM2, Cov_BLM2, Wid_BLM2, Save_cce, Cov_cce, Wid_cce, Save_Bai13, Cov_Bai13, Wid_Bai13, Save_Bai14, Cov_Bai14, Wid_Bai14)
    }

    # Store Results
    idx <- (n - 1) * length(Tlist) + t
    save_result[idx, ] <- c(
      N,
      T,
      mean(res[, 1]) - 1,                      # Bias
      mean((res[, 1] - mean(res[, 1]))^2),     # Variance
      mean(res[, 2]),                          # Coverage
      mean(res[, 3]),                          # Width
      mean(res[, 4]),
      mean(res[, 5]),
      mean(res[, 6]) - 1,                      # Bias CF
      mean((res[, 6] - mean(res[, 6]))^2),     # Variance CF
      mean(res[, 7]),                          # Coverage CF
      mean(res[, 8]),                          # Width CF
      mean(res[, 9]) - 1,                      # Bias TWFE
      mean((res[, 9] - mean(res[, 9]))^2),     # Variance TWFE
      mean(res[, 10]),                         # Coverage TWFE
      mean(res[, 11]),                         # Width TWFE
      mean(res[, 12]) - 1,                     # Bias Bai
      mean((res[, 12] - mean(res[, 12]))^2),   # Variance Bai
      mean(res[, 13]),                         # Coverage Bai
      mean(res[, 14]),                         # Width Bai
      mean(res[, 15]) - 1,                     # Bias GFE
      mean((res[, 15] - mean(res[, 15]))^2),   # Variance GFE
      mean(res[, 16]),                         # Coverage GFE
      mean(res[, 17]),                         # Width GFE
      mean(res[, 18]) - 1,                     # Bias FA
      mean((res[, 18] - mean(res[, 18]))^2),   # Variance FA
      mean(res[, 19]),                         # Coverage FA
      mean(res[, 20]),                         # Width FA
      mean(res[, 21]) - 1,                     # Bias BLM1
      mean((res[, 21] - mean(res[, 21]))^2),   # Variance BLM1
      mean(res[, 22]),                         # Coverage BLM1
      mean(res[, 23]),                         # Width BLM1
      mean(res[, 24]) - 1,                     # Bias BLM2
      mean((res[, 24] - mean(res[, 24]))^2),   # Variance BLM2
      mean(res[, 25]),                         # Coverage BLM2
      mean(res[, 26]),                         # Width BLM2
      mean(res[, 27]) - 1,                     # Bias CCE
      mean((res[, 27] - mean(res[, 27]))^2),   # Variance CCE
      mean(res[, 28]),                         # Coverage CCE
      mean(res[, 29]),                         # Width CCE
      mean(res[, 30]) - 1,                     # Bias Bai13
      mean((res[, 30] - mean(res[, 30]))^2),   # Variance Bai13
      mean(res[, 31]),                         # Coverage Bai13
      mean(res[, 32]),                         # Width Bai13
      mean(res[, 33]) - 1,                     # Bias Bai14
      mean((res[, 33] - mean(res[, 33]))^2),   # Variance Bai14
      mean(res[, 34]),                         # Coverage Bai14
      mean(res[, 35])                          # Width Bai14
    )
  }
}

# Convert results to data frame
results <- data.frame(save_result)

# Set column names
colnames(results) <- c(
  "N", "T", "Bias", "Var", "Cov", "Wid", "G", "C",
  "Bias CF", "Var CF", "Cov CF", "Wid CF",
  "Bias TWFE", "Var TWFE", "Cov TWFE", "Wid TWFE",
  "Bias Bai", "Var Bai", "Cov Bai", "Wid Bai",
  "Bias GFE", "Var GFE", "Cov GFE", "Wid GFE",
  "Bias FA", "Var FA", "Cov FA", "Wid FA",
  "Bias BLM1", "Var BLM1", "Cov BLM1", "Wid BLM1",
  "Bias BLM2", "Var BLM2", "Cov BLM2", "Wid BLM2",
  "Bias CCE", "Var CCE", "Cov CCE", "Wid CCE",
  "Bias Bai13", "Var Bai13", "Cov Bai13", "Wid Bai13",
  "Bias Bai14", "Var Bai14", "Cov Bai14", "Wid Bai14"
)

# End Timer
toc()

#filename <- paste("results", "_DGP" , DGP, "_rho",rho,"_kappa", kappa ,".RData",sep="")

#save(results, file=filename)


results

toc()
