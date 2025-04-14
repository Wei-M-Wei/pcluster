#' Estimator Function with Optional Cross-Fitting
#'@title Inference after discretizing unobserved heterogeneity
#'
#' @description R package 'pcluster' is dedicated to do inference for linear panel data model with unkown functions of fixed effects.
#'
#' @param formula regression formula
#' @param data data which contains id, time, Y, X
#' @param index name of 'id' and 'time' in the data
#' @param CF cross-fitting if true
#' @param init initial iteration number of k-means
#'
#' @returns A list of fitted results is returned.
#' Within this outputted list, the following elements can be found:
#'     \item{res}{regression model.}
#'     \item{G}{cluster number of unit group.}
#'     \item{C}{cluster number of time group.}
#'     \item{estimate_correct}{corrected standard error, t value, and p value.}
#'     \item{summary_table}{summary of the model together with corrected estimates in the 'Coefficients'.}
#'
#' @import plm
#' @export
#'
#' @examples set.seed(1)

#' # Generate data

#' N <- 50
#' T <- 50
#' E <- matrix( rnorm(N*T), N, T)
#' U <- matrix( rnorm(N*T), N, T)
#' E[, 1] <- rnorm(N)
#' U[, 1] <- rnorm(N)
#'
#'
#' A <- rgamma(N, 1)
#' B <- rgamma(N, 1)
#'
#' A_rep <- matrix(rep(A, T), N, T)
#' B_rep <- t(matrix(rep(B, N), T, N))
#'
#'
#' F <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 10)
#' H <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 5)
#'
#'
#' X <- H + U
#' Y <- X + F + E
#'
#' id_code <- numeric(N * T)
#' for (i in 1:N) {
#'   for (t in 1:T) {
#'     id_code[((t - 1) * N + i):(t * N)] <- 10+i
#'   }
#' }
#'
#' time_code <- numeric(N * T)
#' for (i in 1:N) {
#'   for (t in 1:T) {
#'     time_code[((t - 1) * N + 1):(t * N)] <- 10 + t
#'   }
#' }
#'
#' vY <- c(Y)
#' vX <- c(X)
#'
#' data <- data.frame(id_code, time_code, vY, vX)
#' formula <- vY ~ vX - 1
#' index = c("id_code", "time_code")
#' init <- 300
#'
#' # Baseline estimate
#' est <- estimator_dc(formula, data, index, init = init)
#' ols <- est[["res"]]
#' G <- est[["G"]]
#' C <- est[["C"]]
#'
#' # 'summary_correct' contains the original estimates and standard error corrected by the heteroskedasticity autocorrelation consistent standard error
#' summary_correct = est$summary_table
#' coef_estimate = summary_correct$coefficients$Estimate
#' std_error_original = summary_correct$coefficients$`Std. Error`
#' std_error_corrected = summary_correct$coefficients$`Std. Error corrected`
#'
#' # Cross-fitted estimate
#' est_CF <- estimator_dc(formula, data, index, CF = TRUE, init = init)
#' ols_CF <- est_CF[["res"]]
#' summary_correct_CF = est_CF$summary_table
#' coef_estimate_CF = summary_correct_CF$coefficients$Estimate
#' std_error_original_CF = summary_correct_CF$coefficients$`Std. Error`
#' std_error_corrected_CF = summary_correct_CF$coefficients$`Std. Error corrected`

estimator_dc <- function(formula, data, index, CF = FALSE, init = 30) {
  if (!CF) {
    est_without_CF(formula, data, index, init)

  } else {
    est_with_CF(formula, data,index, init)
  }
}

#' General Clustering Function
#' @param Y outcome variable
#' @param X_list covariate matrix
#' @param N sample size
#' @param T panel size
#' @param init initial iteration number of k-means
#' @param type cluster across units if 'tall' or time periods if 'long'
#' @param groups number of the cluster and the default is 'NULL'
#'
#' @returns A list of fitted results is returned.
#' Within this outputted list, the following elements can be found:
#'     \item{res}{cluster results.}
#'     \item{clusters}{cluster number.}
#'
#' @export
cluster_general <- function(Y, X_list, N, T, init, type = "long", groups = NULL) {
  if (type == "long") {
    Y_mean <- rowMeans(Y)
    X_means <- lapply(X_list, rowMeans)
    data <- cbind(Y_mean, do.call(cbind, X_means))
    variance <- sum(sapply(1:T, function(t) {
      norm(cbind(Y[, t], sapply(X_list, function(X) X[, t])) - data, type = "F")^2
    })) / (N * T^2)
    dim_size <- N
  } else if (type == "tall") {
    Y_mean <- colMeans(Y)
    X_means <- lapply(X_list, colMeans)
    data <- cbind(Y_mean, do.call(cbind, X_means))
    variance <- sum(sapply(1:N, function(i) {
      norm(cbind(Y[i, ], sapply(X_list, function(X) X[i, ])) - data, type = "F")^2
    })) / (N^2 * T)
    dim_size <- T
  } else {
    stop("Invalid type. Use 'long' for rows or 'tall' for columns.")
  }

  # Clustering logic
  if (!is.null(groups)) {
    clusters <- groups
    k_result <- kmeans(data, centers = clusters, algorithm = "Hartigan-Wong", nstart = init)
  } else {
    clusters <- 1
    repeat {
      k_result <- kmeans(data, centers = clusters, algorithm = "Hartigan-Wong", nstart = init)
      if (k_result$tot.withinss / dim_size <= variance) break
      clusters <- clusters + 1
    }
  }
  list(res = k_result, clusters = clusters)
}

# Helper Function to Reshape Data into a Matrix
reshape_to_matrix <- function(data, id_var, time_var, value_var) {
  reshaped <- reshape(data[, c(id_var, time_var, value_var)],
                      idvar = id_var,
                      timevar = time_var,
                      direction = "wide")
  as.matrix(reshaped[, -1])  # Drop the 'id' column and return the matrix
}


# Function to recode indices
recode_indices <- function(data, index) {

  unique_id <- sort(unique(data[,c(index[1])]))
  index_to_id <- setNames(1:length(unique_id), unique_id)
  data$id <- index_to_id[as.character(data[,c(index[1])])]

  unique_time <- sort(unique(data[,c(index[2])]))
  index_to_time <- setNames(1:length(unique_time), unique_time)
  data$time <- index_to_time[as.character(data[,c(index[2])])]

  data
}

# Function to Split Data into Dynamic Folds
split_indices <- function(dim, folds) {
  indices <- seq_len(dim)
  split(seq_len(dim), cut(indices, folds, labels = FALSE))
}

# Estimator Without Cross-Fitting
est_without_CF <- function(formula, data, index, init) {
  formula_vars <- all.vars(formula)
  dependent_var <- formula_vars[1]
  independent_vars <- formula_vars[-1]

  data <- recode_indices(data,index)

  Y <- reshape_to_matrix(data, "id", "time", dependent_var)
  X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
  X_combined <- do.call(cbind, X_list)

  N <- nrow(Y)
  T <- ncol(Y)

  clusteri <- cluster_general(Y, X_list, N, T, init, type = "long")
  G <- clusteri$clusters
  klong <- clusteri$res


  clustert <- cluster_general(Y, X_list, N, T, init, type = "tall")
  C <- clustert$clusters
  ktall <- clustert$res

  Du <- matrix(0, N, G)
  Dv <- matrix(0, T, C)

  for (j in seq_len(G)) {
    Du[, j] <- as.numeric(klong$cluster == j)
  }
  for (j in seq_len(C)) {
    Dv[, j] <- as.numeric(ktall$cluster == j)
  }



  Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
  Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)

  tY <- c(Mu %*% Y %*% Mv)
  tX_combined <- do.call(cbind, lapply(X_list, function(X) c(Mu %*% X %*% Mv)))
  new_data <- data.frame(tY, tX_combined,  data[, "id"], data[, "time"])
  colnames(new_data) <- c(formula_vars, "id", "time")
  res = plm(formula, data = new_data, index = c("id", "time"), model = "pooling")

  coefs <- coef(res)
  se_corrected <- sqrt(vcovHC(res, type = "HC0", method = "arellano")) * sqrt((N * T) / ((N-G)*(T-C)))
  t_values_corrected <- coefs / se_corrected

  # Calculate p-values from t-distribution for each coefficient
  df <- res$df.residual  # degrees of freedom
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    SE_Corrected = se_corrected,
    t_value_Corrected = t_values_corrected,
    p_value_Corrected = p_values_corrected
  )

  colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

  summary_table_correct$Signif <- sapply(summary_table_correct[["Pr(>|t|) corrected"]], get_stars)

  summary_table <- list(
    call = summary(res)$call,
    coefficients = summary_table_correct
  )
  summary_table$significance_codes <- "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

  list(res = res, G = G, C = C, estimate_correct = summary_table_correct, summary_table = summary_table)
}

# Estimator with Cross-Fitting
est_with_CF <- function(formula, data, index, init, folds = 2) {
  formula_vars <- all.vars(formula)
  dependent_var <- formula_vars[1]
  independent_vars <- formula_vars[-1]

  data <- recode_indices(data,index)

  Y <- reshape_to_matrix(data, "id", "time", dependent_var)
  X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))

  N <- nrow(Y)
  T <- ncol(Y)
  N_folds <- split_indices(N, folds)
  T_folds <- split_indices(T, folds)

  fold_combinations <- expand.grid(N_fold = seq_along(N_folds), T_fold = seq_along(T_folds))

  Gd <- numeric(nrow(fold_combinations))
  Cd <- numeric(nrow(fold_combinations))
  projection_matrices <- list()

  for (i in seq_len(nrow(fold_combinations))) {
    N_idx <- N_folds[[fold_combinations$N_fold[i]]]
    T_idx <- T_folds[[fold_combinations$T_fold[i]]]

    dim_N <- length(N_idx)
    dim_T <- length(T_idx)

    Y_comp_long <- Y[N_idx, -T_idx]
    X_comp_long_list <- lapply(X_list, function(X) X[N_idx, -T_idx])

    Y_comp_tall <- Y[-N_idx, T_idx]
    X_comp_tall_list <- lapply(X_list, function(X) X[-N_idx, T_idx])

    clusteri <- cluster_general(Y_comp_long, X_comp_long_list, nrow(Y_comp_long), ncol(Y_comp_long), init, type = "long")
    clustert <- cluster_general(Y_comp_tall, X_comp_tall_list, nrow(Y_comp_tall), ncol(Y_comp_tall), init, type = "tall")

    G <- clusteri$clusters
    C <- clustert$clusters

    Du <- matrix(0, dim_N, G)
    Dv <- matrix(0, dim_T, C)

    for (j in seq_len(G)) {
      Du[, j] <- as.numeric(clusteri$res$cluster == j)
    }
    for (j in seq_len(C)) {
      Dv[, j] <- as.numeric(clustert$res$cluster == j)
    }

    Mu <- diag(1, dim_N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
    Mv <- diag(1, dim_T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)

    projection_matrices[[i]] <- list(Mu = Mu, Mv = Mv, N_slice = N_idx, T_slice = T_idx)
    Gd[i] <- G
    Cd[i] <- C
  }

  df <- 0

  projection <- function(X) {
    for (i in seq_along(projection_matrices)) {
      proj <- projection_matrices[[i]]
      X[proj$N_slice, proj$T_slice] <- proj$Mu %*% X[proj$N_slice, proj$T_slice] %*% proj$Mv
    }
    X
  }

  Y <- projection(Y)
  X_list <- lapply(X_list, projection)

  for (i in seq_along(projection_matrices)) {
    proj <- projection_matrices[[i]]
    df <- df + (length(proj$N_slice) - Gd[i]) * (length(proj$T_slice) - Cd[i])
  }

  tY <- c(Y)
  tX_combined <- do.call(cbind, lapply(X_list, function(X) c(X)))
  new_data <- data.frame(tY, tX_combined, data[, "id"], data[, "time"])
  colnames(new_data) <- c(formula_vars, "id", "time")

  res = plm(formula, data = new_data, index = c("id", "time"), model = "pooling")

  coefs <- coef(res)
  se_corrected <- sqrt(vcovHC(res, type = "HC0", method = "arellano")) * sqrt((N * T) / df)
  t_values_corrected <- coefs / se_corrected

  # Calculate p-values from t-distribution for each coefficient
  df_res <- res$df.residual  # degrees of freedom
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df_res)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    SE_Corrected = se_corrected,
    t_value_Corrected = t_values_corrected,
    p_value_Corrected = p_values_corrected
  )

  colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

  summary_table_correct$Signif <- sapply(summary_table_correct[["Pr(>|t|) corrected"]], get_stars)

  summary_table <- list(
    call = summary(res)$call,
    coefficients = summary_table_correct,
    significance_codes = "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
  )


  list(res = res, df = df, estimate_correct = summary_table_correct, summary_table = summary_table)
}

transform_vector <- function(vec, N) {
  C <- length(vec)
  mat <- matrix(0, nrow = N, ncol = N * C)

  for (i in 1:N) {
    start_pos <- (i - 1) * C + 1
    end_pos <- i * C
    mat[i, start_pos:end_pos] <- vec
  }

  return(mat)
}

transform_matrix <- function(mat, T, t) {
  N <- nrow(mat)
  G <- ncol(mat)

  # Initialize an N x (T * G) matrix with zeros
  new_mat <- matrix(0, nrow = N, ncol = T * G)

  # Compute the column range for insertion
  start_col <- (t - 1) * G + 1
  end_col <- t * G

  # Insert the original matrix into the computed column range
  new_mat[, start_col:end_col] <- mat

  return(new_mat)
}

get_stars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else if (p < 0.1) return(".")
  else return(" ")
}


# est_without_CF_dummy <- function(formula, data, index, init) {
#   formula_vars <- all.vars(formula)
#   dependent_var <- formula_vars[1]
#   independent_vars <- formula_vars[-1]
#
#   data <- recode_indices(data,index)
#
#   Y <- reshape_to_matrix(data, "id", "time", dependent_var)
#   X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
#   X_combined <- do.call(cbind, X_list)
#
#   N <- nrow(Y)
#   T <- ncol(Y)
#
#   clusteri <- cluster_general(Y, X_list, N, T, init, type = "long")
#   G <- clusteri$clusters
#   klong <- clusteri$res
#
#   clustert <- cluster_general(Y, X_list, N, T, init, type = "tall")
#   C <- clustert$clusters
#   ktall <- clustert$res
#
#   Du <- matrix(0, N, G)
#   Dv <- matrix(0, T, C)
#
#   for (j in seq_len(G)) {
#     Du[, j] <- as.numeric(klong$cluster == j)
#   }
#   for (j in seq_len(C)) {
#     Dv[, j] <- as.numeric(ktall$cluster == j)
#   }
#
#   sig_dummy = matrix(0,N*T,N*C)
#   for (t in 1:T) {
#     sig_dummy[((t-1)*N+1):(t*N), ] = transform_vector(Dv[t,],N)
#
#   }
#   mu_dummy = matrix(0,N*T,T*G)
#   for (t in 1:T) {
#     mu_dummy[((t-1)*N+1):(t*N), ] = transform_matrix(Du,T,t)
#
#   }
#
#   new_data <- data.frame(data$vY, data$vX, sig_dummy, mu_dummy,   data[, "id"], data[, "time"])
#   vec = colnames(new_data)[-1]
#   formula_str <- paste(colnames(new_data)[1], "~-1+", paste(vec[-((length(vec)-1):length(vec))], collapse = " + "))
#
#   colnames(new_data)[length(colnames(new_data))-1] <- c( 'id')
#   colnames(new_data)[length(colnames(new_data))] <- c( 'time')
#   res_dummy = plm(formula_str, data = new_data, index = c("id", "time"),  model = "pooling")
#   list(res = res_dummy, G = G, C = C)
# }
#
# estimator_dc_dummy <- function(formula, data, index, CF = FALSE, init = 30) {
#   if (!CF) {
#     est_without_CF_dummy(formula, data, index, init)
#   }
# }
