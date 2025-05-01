# Estimator Function with Clustering
estimator_dci <- function(formula, data) {
  # Extract variable names from the formula
  formula_vars <- all.vars(formula)
  dependent_var <- formula_vars[1]
  independent_vars <- formula_vars[-1]

  # Ensure the required columns exist in the data
  if (!all(c("id", "time", dependent_var, independent_vars) %in% colnames(data))) {
    stop("Data must contain 'id', 'time', dependent, and independent variables.")
  }

  # Perform clustering for rows and columns
  N <- length(unique(data$id))
  T <- length(unique(data$time))
  K <- length(formula_vars) - 1
  hat_beta <- matrix(0, N, K)

  for (i in seq_len(N)) {
    Y <- reshape_to_matrix(data[data$id != i, ], "id", "time", dependent_var)
    X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data[data$id != i, ], "id", "time", var))
    clustert <- cluster_general(Y, X_list, N - 1, T, type = "tall")
    C <- clustert$clusters
    ktall <- clustert$res
    Dv <- matrix(0, T, C)

    for (j in seq_len(C)) {
      Dv[, j] <- as.numeric(ktall$cluster == j)
    }

    Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
    tY <- Mv %*% as.vector(data[data$id == i, dependent_var])
    X_list_i <- lapply(independent_vars, function(var) reshape_to_matrix(data[data$id == i, ], "id", "time", var))
    tX_combined <- do.call(cbind, lapply(X_list_i, function(X) Mv %*% as.vector(X)))
    new_data <- data.frame(tY, tX_combined, data[, "id"], data[, "time"])
    colnames(new_data) <- c(formula_vars, "id", "time")
    hat_beta[i, ] <- lm(formula, data = new_data)$coefficients
  }

  # Perform OLS regression and return results
  list(beta = mean(hat_beta), beta_var = var(hat_beta))
}


# Helper Function to Reshape Data into a Matrix
reshape_to_matrix <- function(data, id_var, time_var, value_var) {
  reshaped <- reshape(data[, c(id_var, time_var, value_var)],
                      idvar = id_var,
                      timevar = time_var,
                      direction = "wide")
  as.matrix(reshaped[, -1])  # Drop the 'id' column and return the matrix
}

# Data Generation Function
generatei <- function(N, T, DGP, rho = 0, kappa = 0) {
  # Generate Random Error Matrices and Unobservables
  E <- matrix(0, N, T)
  U <- matrix(0, N, T)
  E[, 1] <- rnorm(N)
  U[, 1] <- rnorm(N)

  for (t in 2:T) {
    E[, t] <- kappa * E[, t - 1] + sqrt(1 - kappa^2) * rnorm(N)
    U[, t] <- kappa * U[, t - 1] + sqrt(1 - kappa^2) * rnorm(N)
  }

  A <- rgamma(N, 1)
  BURN <- 10000
  B0 <- rgamma(1, 1)
  for (t in 1:BURN) {
    B1 <- rho * B0 + rgamma(1, (1 - rho) * (1 - rho) / (1 - rho^2), rate = (1 - rho) / (1 - rho^2))
    B0 <- B1
  }
  B <- numeric(T)
  B[1] <- B0
  for (t in 2:T) {
    B[t] <- rho * B[t - 1] + rgamma(1, (1 - rho) * (1 - rho) / (1 - rho^2), rate = (1 - rho) / (1 - rho^2))
  }

  # Compute H and G Matrices
  A_rep <- matrix(rep(A, T), N, T)
  B_rep <- t(matrix(rep(B, N), T, N))
  if (DGP == 1) {
    F <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 10)
    H <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 5)
  } else if (DGP == 2) {
    F <- A_rep^2 + B_rep * A_rep + sin(B_rep * A_rep)
    H <- B_rep^2 + B_rep * A_rep + sin(B_rep * A_rep)
  }

  # Generate X and Y Matrices
  beta <- 4 * runif(N)
  X <- H + U
  Y <- matrix(0, N, T)
  for (i in 1:N) {
    Y[i, ] <- beta[i] * X[i, ] + F[i, ] + E[i, ]
  }

  id <- numeric(N * T)
  for (i in 1:N) {
    for (t in 1:T) {
      id[(t - 1) * N + i] <- i
    }
  }

  time <- numeric(N * T)
  for (i in 1:N) {
    for (t in 1:T) {
      time[((t - 1) * N + 1):(t * N)] <- t
    }
  }
  vY <- c(Y)
  vX <- c(X)

  list(data = data.frame(id, time, vY, vX), beta = beta)
}


