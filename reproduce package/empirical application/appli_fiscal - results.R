# Clear workspace
rm(list = ls())

# Source external functions
source("functions.R")

# Load required libraries
library(readr)
library(sandwich)
library(gtools)
library(plm)
library(pcluster)
library(dplyr)
library(kableExtra)

# Load and preprocess data
df <- read_csv("57.to.2008.raw.data.Aug4.2014.csv") %>%
  filter(year > 1957, state != 2) %>%
  arrange(year, state)

# Initialize matrices to store results
save_result <- matrix(0, nrow = 10, ncol = 6)
save_result_Groups <- matrix(0, nrow = 2, ncol = 5)
pvalues <- matrix(0, nrow = 5, ncol = 5)

# Define constants
N <- 49
T <- 51

# Loop over different outcomes
for (ind in 1:5) {
  # Select dependent variable based on index
  Y <- switch(ind,
              `1` = df$`inc tax pi`,
              `2` = df$nonresrevpi,
              `3` = df$totexppi,
              `4` = df$`edu exp pi`,
              `5` = df$`savings pi`)

  X <- df$resrevpi

  # Create ID and time vectors
  id <- numeric(N * T)
  for (i in 1:N) {
    for (t in 1:T) {
      id[(t - 1) * N + i] <- i
    }
  }

  time <- numeric(N*T)
  for (i in (1:N)){
    for (t in (1:T)){
      time[((t-1)*N+1):(t*N)]<-t
    }
  }

  # Prepare data frame
  data <- data.frame(
    Y = c(as.matrix(Y)),
    X = c(as.matrix(X)),
    id = id,
    time = time
  )

  # TWFE estimate
  ols_twfe <- plm(Y ~ X, data = data, index = c("id", "time"), model = "within", effect = "twoways")
  std_errors_twfe <- sqrt(vcovHC(ols_twfe, type = "HC0", method = "arellano"))
  Save_twfe <- ols_twfe$coefficients

  # Standardize Y and X
  sY <- sd(data$Y)
  sX <- sd(data$X)
  data <- data %>%
    mutate(
      Y = (Y - mean(Y)) / sY,
      X = (X - mean(X)) / sX
    )

  # Baseline estimate
  set.seed(1)
  est <- estimator_dc(Y ~ X - 1, data, index = c('id', 'time'))
  ols <- est[["res"]]
  G <- est[["G"]]
  C <- est[["C"]]
  std_errors <- est$summary_table$coefficients$`Std. Error corrected` * sY / sX
  Save <- ols$coefficients * sY / sX

  # Bai's estimate with R = T^{1/4}
  mY <- reshape_to_matrix(data, "id", "time", "Y")
  X_list <- lapply(c("X"), function(var) reshape_to_matrix(data, "id", "time", var))
  mX <- do.call(cbind, X_list)
  R <- floor(T^(1/4))
  est_Bai <- LS.factor(mX, mY, R = R)
  Save_Bai <- est_Bai[["beta"]] * sY / sX
  std_errors_Bai <- est_Bai[["hac.se"]] * sY / sX

  # FA estimate
  est_FA <- FA(Y ~ X - 1, data)
  ols_FA <- est_FA[["res"]]
  std_errors_FA <- sqrt(vcovHC(ols_FA, type = "HC0", method = "arellano")) * sY / sX
  Save_FA <- ols_FA$coefficients * sY / sX

  # CCE estimator
  pdata <- pdata.frame(data, index = c("id", "time"))
  ccep <- pcce(Y ~ X, data = pdata, model = "p")
  std_errors_cce <- sqrt(vcovHC(ccep, type = "HC0", method = "arellano")) * sY / sX
  Save_cce <- ccep$coefficients * sY / sX

  # Save results
  save_result[2 * (ind - 1) + 1, 2:6] <- c(Save, Save_twfe, Save_Bai, Save_FA, Save_cce)
  save_result[2 * ind, 2:6] <- c(std_errors, std_errors_twfe, std_errors_Bai, std_errors_FA, std_errors_cce)
  save_result_Groups[, ind] <- c(G, C)
  pvalues[ind, ] <- (1 - c(pnorm(abs(Save) / std_errors), pnorm(abs(Save_twfe) / std_errors_twfe),
                           pnorm(abs(Save_Bai) / std_errors_Bai), pnorm(abs(Save_FA) / std_errors_FA),
                           pnorm(abs(Save_cce) / std_errors_cce)))*2
}

# Prepare results for output
save_result <- as.data.frame(save_result)

colnames(save_result_Groups) <- c("Nonresource revenue", "Income tax revenue", "Total expenditures",
                                  "Education expenditures", "Savings")
rownames(save_result_Groups) <- c("G", "C")
save_result_Groups <- as.data.frame(save_result_Groups)

# Round and output results as LaTeX tables
save_result %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  mutate(across(where(is.numeric), ~ ifelse(row_number() %% 2 == 0, paste0("(", ., ")"), .))) %>%
  kable("latex", align = "c", linesep = "", escape = FALSE, row.names = FALSE)

save_result_Groups %>%
  kable("latex", align = "c", linesep = "", escape = FALSE)

# Function to add significance stars
add_stars <- function(p) {
  ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
}

# Apply the function to the pvalues data frame and create a data frame of stars
stars <- apply(pvalues, 1:2, add_stars)

# Combine the save_result with stars, pasting stars to the corresponding values
save_result_with_stars <- save_result

# Round the results and prepare for LaTeX output
save_result_with_stars <- save_result_with_stars %>%
  mutate(across(everything(), ~ round(as.numeric(gsub("[*]+", "", .)), 3))) %>%
  mutate(across(everything(), ~ ifelse(row_number() %% 2 == 0, paste0("(", ., ")"), .)))

save_result_with_stars[1, -1] <- paste0(save_result_with_stars[1, -1], stars[1, ])
save_result_with_stars[3, -1] <- paste0(save_result_with_stars[3, -1], stars[2, ])
save_result_with_stars[5, -1] <- paste0(save_result_with_stars[5, -1], stars[3, ])
save_result_with_stars[7, -1] <- paste0(save_result_with_stars[7, -1], stars[4, ])
save_result_with_stars[9, -1] <- paste0(save_result_with_stars[9, -1], stars[5, ])

save_result_with_stars[, 1] <- c("Nonresource revenue", " ", "Income tax revenue", " ", "Total expenditures", " ",
                                 "Education expenditures", " ", "Savings", " ")

# Output the modified results as a LaTeX table
save_result_with_stars %>%
  kable("latex", align = "c", linesep = "", escape = FALSE)
