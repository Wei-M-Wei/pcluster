# Helper Function to Reshape Data into a Matrix
reshape_to_matrix <- function(data, id_var, time_var, value_var) {
  reshaped <- reshape(data[, c(id_var, time_var, value_var)],
                      idvar = id_var,
                      timevar = time_var,
                      direction = "wide")
  as.matrix(reshaped[, -1])  # Drop the 'id' column and return the matrix
}

# Function to Split Data into Dynamic Folds
split_indices <- function(dim, folds) {
  indices <- seq_len(dim)
  split(seq_len(dim), cut(indices, folds, labels = FALSE))
}


# Data Generation Function
generate <- function(N, T, DGP, rho = 0, kappa = 0) {
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

  A_rep <- matrix(rep(A, T), N, T)
  B_rep <- t(matrix(rep(B, N), T, N))

  if (DGP == 1) {
    F <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 10)
    H <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 5)
  } else if (DGP == 2) {
    F <- A_rep^2 + B_rep * A_rep + sin(B_rep * A_rep)
    H <- B_rep^2 + B_rep * A_rep + sin(B_rep * A_rep)
  }

  X <- H + U
  Y <- X + F + E

  id <- numeric(N * T)
  for (i in 1:N) {
    for (t in 1:T) {
      id[((t - 1) * N + i):(t * N)] <- i
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

  data.frame(id, time, vY, vX)
}

LS.factor<- function(X,Y,R = 3,tol = 0.0001,repMin = 20,repMax = 100,exit_tol = 0.1){
  N = dim(Y)[1]
  T = dim(Y)[2]
  XX = X
  if (is.na(dim(X)[3])){XX = array(X,dim=c(N,T,1))}
  K = dim(XX)[3]


  beta_init = as.matrix(rnorm(n=K), dim = c(K,1))
  Fhat_init = array(0, dim = c(T,R))
  obj_init = 0

  beta_hat_int = beta_init + 1
  Fhat_interim = Fhat_init + 1
  obj_int = obj_init + 1
  iterations = 0

  X.tilde = array(NA, dim = c(N,T, K))

  Y.squeeze<-as.matrix(as.vector(Y),dim = N*T)
  X.squeeze<-matrix(XX, nrow = N*T, ncol = K)

  W.squeeze = Y.squeeze - X.squeeze%*%beta_init
  Lambda_hat_int = array(NA, dim = c(N,R))
  for (rr in 1: R){
    Lambda_hat_int[,rr] = rnorm(n = N)
  }
  Lambda_init = Lambda_hat_int - 1

  ef = 1
  while(
    abs(obj_init-obj_int)>=tol
    &(ef==1|iterations<repMin)
    &iterations<repMax
  ){
    iterations = iterations + 1
    Mf = diag(T) - Fhat_interim%*%t(Fhat_interim)/T
    Mlambda = diag(N) - (Lambda_hat_int)%*%solve(t(Lambda_hat_int)%*%Lambda_hat_int)%*%t(Lambda_hat_int)

    for (k in 1:K){
      X.tilde[,,k] = Mlambda%*%(XX[,,k])%*%Mf
    }

    Y.tilde = Mlambda%*%(Y)%*%Mf

    Y.squeeze.tilde<-as.matrix(as.vector(Y.tilde),dim = N*T)
    X.squeeze.tilde<-matrix(X.tilde, nrow = N*T, ncol = K)
    X.squeeze<-matrix(X, nrow = N*T, ncol = K)

    beta_init = beta_hat_int
    beta_hat_int = solve(t(X.squeeze)%*%X.squeeze.tilde)%*%t(X.squeeze)%*%Y.squeeze.tilde

    W.squeeze = Y.squeeze - X.squeeze%*%beta_init
    W = matrix(W.squeeze,nrow = N,ncol = T)

    Fhat_init = Fhat_interim
    Fhat_interim = sqrt(T)*as.matrix(eigen(t(W)%*%(W))$vectors[,1:R])

    Lambda_init = Lambda_hat_int
    Lambda_hat_int = ((1/T)*W%*%(Fhat_interim))

    obj_init = obj_int
    obj_int <- sum(diag((W - Lambda_hat_int%*%t(Fhat_interim))%*%
                          t(W - Lambda_hat_int%*%t(Fhat_interim))))
    error<- W - Lambda_hat_int%*%t(Fhat_interim)

    svd<-svd(W)

    if(iterations>repMin&obj_int>obj_init*(1+exit_tol)){ef = -1} else{ef = 1}

    if(ef == -1){beta_init = beta_hat_int}
  }

  s.square <- sum((Y.squeeze.tilde - X.squeeze.tilde%*%beta_hat_int)^2)/((N - R)*(T - R) - 1)
  u <- Y.squeeze.tilde - X.squeeze.tilde%*%beta_hat_int
  X.X <- t(X.squeeze.tilde)%*%X.squeeze.tilde
  xs <- apply(X.squeeze.tilde,MARGIN = 2, function(x){u*x})
  beta.hac.se <- diag(sqrt((N*T)/((N - R)*(T - R) - 1)*solve(X.X)%*%t(xs)%*%xs%*%solve(X.X)))
  beta.se     <- diag(sqrt(s.square*solve(X.X)))


  result<- list(beta = beta_hat_int, hac.se = beta.hac.se, se = beta.se,
                F = Fhat_interim, Lambda = Lambda_hat_int,
                error = error,svd = svd, Iterations = iterations, Exit_flag = ef)

  return(result)
}

hnn.est <- function(Y,X, L, F){
  N<-dim(Y)[1]
  T<-dim(Y)[2]
  L <- as.matrix(L)
  F <- as.matrix(F)
  deltalambda<-array(NA, dim = c(N,N))
  for (i in 1:N){
    for (j in 1:N){
      deltalambda[i,j]<-sum((L[i,] - L[j,])^2)
    }
  }
  diag(deltalambda) <- Inf
  # deltalambda[!lower.tri(deltalambda)] <- Inf

  check <- N
  hnn.i <- list(c(1:N))
  hnn.inner <- c(which(deltalambda == min(deltalambda), arr.ind = TRUE)[1,])
  hnn.i[[1]] <- hnn.i[[1]]
  hnn.i[[2]] <- append(hnn.i[[1]][-hnn.inner],list(hnn.inner))
  k <- 2
  while(check != 0){
    singles <- unlist(hnn.i[[k]][which(lengths(hnn.i[[k]]) == 1)])
    twos    <- hnn.i[[k]][which(lengths(hnn.i[[k]]) == 2)]
    threes  <- hnn.i[[k]][which(lengths(hnn.i[[k]]) == 3)]
    fours   <- hnn.i[[k]][which(lengths(hnn.i[[k]]) == 4)]
    if(length(fours)>0){
      four.splits <- fours
      for (kk in 1:length(fours)){
        perms <- gtools::permutations(4,4,fours[[kk]])
        four.inner.dist <- apply(perms, MARGIN = 1, function(x){max(deltalambda[x[1],x[2]],deltalambda[x[3],x[4]])})
        four.splits[[kk]] <- append(list(perms[which.min(four.inner.dist),c(1,2)]),list(perms[which.min(four.inner.dist),c(3,4)]))
      }
      fours <- t(array(unlist(four.splits), dim = c(2,2*length(four.splits))))
      hnn.i[[k]] <- hnn.i[[k]][which(lengths(hnn.i[[k]]) != 4)]
      for (kk in 1:dim(fours)[1]){
        hnn.i[[k]] <- append(hnn.i[[k]], list(fours[kk,]))
      }

    }
    if(length(singles)!=0){
      min.linkage <- as.matrix(deltalambda[,singles])
      hnn.inner   <- c(which(min.linkage == min(min.linkage), arr.ind = TRUE)[1,])
      hnn.inner <- unlist(append(hnn.i[[k]][which(sapply(hnn.i[[k]], function(x){hnn.inner[1] %in% x}))],list(singles[hnn.inner[2]])))
      hnn.remove <- (apply(sapply(hnn.i[[k]], function(x){hnn.inner %in% x}), MARGIN = 2, sum) == 0)
      hnn.i[[k+1]] <- append(subset(hnn.i[[k]], hnn.remove), list(hnn.inner))
    }else{
      hnn.i[[k+1]] <- hnn.i[[k]]
    }
    check <- sum(1 - sapply(hnn.i[[k+1]], function(x){length(x) == 2|length(x) == 3}))
    hnn.i[[k]] <- hnn.i[[k+1]]
  }
  D.N <- array(0, dim = c(N, length(hnn.i[[3]])))
  for (k in 1:length(hnn.i[[3]])){
    D.N[hnn.i[[3]][[k]], k] <- 1
  }
  deltaf<-array(NA, dim = c(T,T))
  for (t in 1:T){
    for (s in 1:T){
      deltaf[t,s]<-sum((F[t,] - F[s,])^2)
    }
  }
  diag(deltaf) <- Inf

  check <- T
  hnn.t <- list(c(1:T))
  hnn.inner <- c(which(deltaf == min(deltaf), arr.ind = TRUE)[1,])
  hnn.t[[1]] <- hnn.t[[1]]
  hnn.t[[2]] <- append(hnn.t[[1]][-hnn.inner],list(hnn.inner))
  k <- 2
  while(check != 0){
    singles <- unlist(hnn.t[[k]][which(lengths(hnn.t[[k]]) == 1)])
    twos    <- hnn.t[[k]][which(lengths(hnn.t[[k]]) == 2)]
    threes  <- hnn.t[[k]][which(lengths(hnn.t[[k]]) == 3)]
    fours   <- hnn.t[[k]][which(lengths(hnn.t[[k]]) == 4)]
    if(length(fours)>0){
      four.splits <- fours
      for (kk in 1:length(fours)){
        perms <- gtools::permutations(4,4,fours[[kk]])
        four.inner.dist <- apply(perms, MARGIN = 1, function(x){max(deltaf[x[1],x[2]],deltaf[x[3],x[4]])})
        four.splits[[kk]] <- append(list(perms[which.min(four.inner.dist),c(1,2)]),list(perms[which.min(four.inner.dist),c(3,4)]))
      }
      fours <- t(array(unlist(four.splits), dim = c(2,2*length(four.splits))))
      hnn.t[[k]] <- hnn.t[[k]][which(lengths(hnn.t[[k]]) != 4)]
      for (kk in 1:dim(fours)[1]){
        hnn.t[[k]] <- append(hnn.t[[k]], list(fours[kk,]))
      }

    }
    if(length(singles)!=0){
      min.linkage <- as.matrix(deltaf[,singles])
      hnn.inner   <- c(which(min.linkage == min(min.linkage), arr.ind = TRUE)[1,])
      hnn.inner <- unlist(append(hnn.t[[k]][which(sapply(hnn.t[[k]], function(x){hnn.inner[1] %in% x}))],list(singles[hnn.inner[2]])))
      hnn.remove <- (apply(sapply(hnn.t[[k]], function(x){hnn.inner %in% x}), MARGIN = 2, sum) == 0)
      hnn.t[[k+1]] <- append(subset(hnn.t[[k]], hnn.remove), list(hnn.inner))
    }else{
      hnn.t[[k+1]] <- hnn.t[[k]]
    }
    check <- sum(1 - sapply(hnn.t[[k+1]], function(x){length(x) == 2|length(x) == 3}))
    hnn.t[[k]] <- hnn.t[[k+1]]
  }
  D.T <- array(0, dim = c(T, length(hnn.t[[3]])))
  for (k in 1:length(hnn.t[[3]])){
    D.T[hnn.t[[3]][[k]], k] <- 1
  }
  t.pair<-apply(D.T, MARGIN = 2, function(x){which(x==1)})

  M.N <- diag(N) - D.N%*%solve(t(D.N)%*%D.N)%*%t(D.N)
  M.T <- diag(T) - D.T%*%solve(t(D.T)%*%D.T)%*%t(D.T)
  Y.tilde = M.N%*%Y%*%M.T
  X.tilde = array(NA, dim = dim(X))
  for (k in 1:dim(X)[3]){
    X.tilde[,,k] = M.N%*%X[,,k]%*%M.T
  }
  Y.squeeze.tilde <- as.matrix(as.vector(Y.tilde),dim = N*T)
  X.squeeze.tilde <- matrix(X.tilde, nrow = N*T, ncol = dim(X)[3])
  X.X <- t(X.squeeze.tilde)%*%X.squeeze.tilde

  beta.hnn <- solve(t(X.squeeze.tilde)%*%X.squeeze.tilde)%*%(t(X.squeeze.tilde)%*%Y.squeeze.tilde)
  # s.square <- sum((Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn)^2)/(N*T - 1)
  se.cluster <- kronecker(D.T,D.N, FUN = '*')
  u <- Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn
  XOX <- array(apply(array(apply(se.cluster,MARGIN = 2, function(x){
    xx <- t(array(X.tilde, dim = c(N*T,dim(X)[3]))[x==1,])%*%diag(u[x==1]^2)%*%array(X.tilde, dim = c(N*T,dim(X)[3]))[x==1,]
  }), dim = c(dim(se.cluster)[2], dim(X)[3])), MARGIN = 2, sum), dim = c(dim(X)[3],dim(X)[3]))
  dfa <- sqrt((N*T - 1)/((N - dim(D.N)[2])*(T - dim(D.T)[2]) - 1))
  beta.cluster.se <- sqrt(solve(X.X)%*%XOX%*%solve(X.X))*dfa

  dfa <- sqrt((N*T - 1)/((N - dim(D.N)[2])*(T - dim(D.T)[2]) - 1))

  s.square <- sum((Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn)^2/(N*T - 1))

  uu <- c(Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn)
  xs <- apply(X.squeeze.tilde,MARGIN = 2, function(x){uu*x})
  beta.hac.se <- sqrt(solve(X.X)%*%(t(xs)%*%xs)%*%solve(X.X))*dfa
  beta.se     <- sqrt(s.square*solve(X.X))*dfa
  return(list(beta.hnn = beta.hnn, hnn.est <- function(Y,X, L, F){
    N<-dim(Y)[1]
    T<-dim(Y)[2]
    L <- as.matrix(L)
    F <- as.matrix(F)
    deltalambda<-array(NA, dim = c(N,N))
    for (i in 1:N){
      for (j in 1:N){
        deltalambda[i,j]<-sum((L[i,] - L[j,])^2)
      }
    }
    diag(deltalambda) <- Inf
    # deltalambda[!lower.tri(deltalambda)] <- Inf

    check <- N
    hnn.i <- list(c(1:N))
    hnn.inner <- c(which(deltalambda == min(deltalambda), arr.ind = TRUE)[1,])
    hnn.i[[1]] <- hnn.i[[1]]
    hnn.i[[2]] <- append(hnn.i[[1]][-hnn.inner],list(hnn.inner))
    k <- 2
    while(check != 0){
      singles <- unlist(hnn.i[[k]][which(lengths(hnn.i[[k]]) == 1)])
      twos    <- hnn.i[[k]][which(lengths(hnn.i[[k]]) == 2)]
      threes  <- hnn.i[[k]][which(lengths(hnn.i[[k]]) == 3)]
      fours   <- hnn.i[[k]][which(lengths(hnn.i[[k]]) == 4)]
      if(length(fours)>0){
        four.splits <- fours
        for (kk in 1:length(fours)){
          perms <- gtools::permutations(4,4,fours[[kk]])
          four.inner.dist <- apply(perms, MARGIN = 1, function(x){max(deltalambda[x[1],x[2]],deltalambda[x[3],x[4]])})
          four.splits[[kk]] <- append(list(perms[which.min(four.inner.dist),c(1,2)]),list(perms[which.min(four.inner.dist),c(3,4)]))
        }
        fours <- t(array(unlist(four.splits), dim = c(2,2*length(four.splits))))
        hnn.i[[k]] <- hnn.i[[k]][which(lengths(hnn.i[[k]]) != 4)]
        for (kk in 1:dim(fours)[1]){
          hnn.i[[k]] <- append(hnn.i[[k]], list(fours[kk,]))
        }

      }
      if(length(singles)!=0){
        min.linkage <- as.matrix(deltalambda[,singles])
        hnn.inner   <- c(which(min.linkage == min(min.linkage), arr.ind = TRUE)[1,])
        hnn.inner <- unlist(append(hnn.i[[k]][which(sapply(hnn.i[[k]], function(x){hnn.inner[1] %in% x}))],list(singles[hnn.inner[2]])))
        hnn.remove <- (apply(sapply(hnn.i[[k]], function(x){hnn.inner %in% x}), MARGIN = 2, sum) == 0)
        hnn.i[[k+1]] <- append(subset(hnn.i[[k]], hnn.remove), list(hnn.inner))
      }else{
        hnn.i[[k+1]] <- hnn.i[[k]]
      }
      check <- sum(1 - sapply(hnn.i[[k+1]], function(x){length(x) == 2|length(x) == 3}))
      hnn.i[[k]] <- hnn.i[[k+1]]
    }
    D.N <- array(0, dim = c(N, length(hnn.i[[3]])))
    for (k in 1:length(hnn.i[[3]])){
      D.N[hnn.i[[3]][[k]], k] <- 1
    }
    deltaf<-array(NA, dim = c(T,T))
    for (t in 1:T){
      for (s in 1:T){
        deltaf[t,s]<-sum((F[t,] - F[s,])^2)
      }
    }
    diag(deltaf) <- Inf

    check <- T
    hnn.t <- list(c(1:T))
    hnn.inner <- c(which(deltaf == min(deltaf), arr.ind = TRUE)[1,])
    hnn.t[[1]] <- hnn.t[[1]]
    hnn.t[[2]] <- append(hnn.t[[1]][-hnn.inner],list(hnn.inner))
    k <- 2
    while(check != 0){
      singles <- unlist(hnn.t[[k]][which(lengths(hnn.t[[k]]) == 1)])
      twos    <- hnn.t[[k]][which(lengths(hnn.t[[k]]) == 2)]
      threes  <- hnn.t[[k]][which(lengths(hnn.t[[k]]) == 3)]
      fours   <- hnn.t[[k]][which(lengths(hnn.t[[k]]) == 4)]
      if(length(fours)>0){
        four.splits <- fours
        for (kk in 1:length(fours)){
          perms <- gtools::permutations(4,4,fours[[kk]])
          four.inner.dist <- apply(perms, MARGIN = 1, function(x){max(deltaf[x[1],x[2]],deltaf[x[3],x[4]])})
          four.splits[[kk]] <- append(list(perms[which.min(four.inner.dist),c(1,2)]),list(perms[which.min(four.inner.dist),c(3,4)]))
        }
        fours <- t(array(unlist(four.splits), dim = c(2,2*length(four.splits))))
        hnn.t[[k]] <- hnn.t[[k]][which(lengths(hnn.t[[k]]) != 4)]
        for (kk in 1:dim(fours)[1]){
          hnn.t[[k]] <- append(hnn.t[[k]], list(fours[kk,]))
        }

      }
      if(length(singles)!=0){
        min.linkage <- as.matrix(deltaf[,singles])
        hnn.inner   <- c(which(min.linkage == min(min.linkage), arr.ind = TRUE)[1,])
        hnn.inner <- unlist(append(hnn.t[[k]][which(sapply(hnn.t[[k]], function(x){hnn.inner[1] %in% x}))],list(singles[hnn.inner[2]])))
        hnn.remove <- (apply(sapply(hnn.t[[k]], function(x){hnn.inner %in% x}), MARGIN = 2, sum) == 0)
        hnn.t[[k+1]] <- append(subset(hnn.t[[k]], hnn.remove), list(hnn.inner))
      }else{
        hnn.t[[k+1]] <- hnn.t[[k]]
      }
      check <- sum(1 - sapply(hnn.t[[k+1]], function(x){length(x) == 2|length(x) == 3}))
      hnn.t[[k]] <- hnn.t[[k+1]]
    }
    D.T <- array(0, dim = c(T, length(hnn.t[[3]])))
    for (k in 1:length(hnn.t[[3]])){
      D.T[hnn.t[[3]][[k]], k] <- 1
    }
    t.pair<-apply(D.T, MARGIN = 2, function(x){which(x==1)})

    M.N <- diag(N) - D.N%*%solve(t(D.N)%*%D.N)%*%t(D.N)
    M.T <- diag(T) - D.T%*%solve(t(D.T)%*%D.T)%*%t(D.T)
    Y.tilde = M.N%*%Y%*%M.T
    X.tilde = array(NA, dim = dim(X))
    for (k in 1:dim(X)[3]){
      X.tilde[,,k] = M.N%*%X[,,k]%*%M.T
    }
    Y.squeeze.tilde <- as.matrix(as.vector(Y.tilde),dim = N*T)
    X.squeeze.tilde <- matrix(X.tilde, nrow = N*T, ncol = dim(X)[3])
    X.X <- t(X.squeeze.tilde)%*%X.squeeze.tilde

    beta.hnn <- solve(t(X.squeeze.tilde)%*%X.squeeze.tilde)%*%(t(X.squeeze.tilde)%*%Y.squeeze.tilde)
    # s.square <- sum((Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn)^2)/(N*T - 1)
    se.cluster <- kronecker(D.T,D.N, FUN = '*')
    u <- Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn
    XOX <- array(apply(array(apply(se.cluster,MARGIN = 2, function(x){
      xx <- t(array(X.tilde, dim = c(N*T,dim(X)[3]))[x==1,])%*%diag(u[x==1]^2)%*%array(X.tilde, dim = c(N*T,dim(X)[3]))[x==1,]
    }), dim = c(dim(se.cluster)[2], dim(X)[3])), MARGIN = 2, sum), dim = c(dim(X)[3],dim(X)[3]))
    dfa <- sqrt((N*T - 1)/((N - dim(D.N)[2])*(T - dim(D.T)[2]) - 1))
    beta.cluster.se <- sqrt(solve(X.X)%*%XOX%*%solve(X.X))*dfa

    dfa <- sqrt((N*T - 1)/((N - dim(D.N)[2])*(T - dim(D.T)[2]) - 1))

    s.square <- sum((Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn)^2/(N*T - 1))

    uu <- c(Y.squeeze.tilde - X.squeeze.tilde%*%beta.hnn)
    xs <- apply(X.squeeze.tilde,MARGIN = 2, function(x){uu*x})
    beta.hac.se <- sqrt(solve(X.X)%*%(t(xs)%*%xs)%*%solve(X.X))*dfa
    beta.se     <- sqrt(s.square*solve(X.X))*dfa
    return(list(beta.hnn = beta.hnn, hac.se = beta.hac.se, se = beta.se, cluster.se = beta.cluster.se,
                D.N = D.N, D.T = D.T,
                Y.diff = Y.squeeze.tilde, X.diff = X.squeeze.tilde))
  } = beta.hac.se, se = beta.se, cluster.se = beta.cluster.se,
  D.N = D.N, D.T = D.T,
  Y.diff = Y.squeeze.tilde, X.diff = X.squeeze.tilde))
}

# Factor Analysis Function
FA <- function(formula, data) {
  formula_vars <- all.vars(formula)
  dependent_var <- formula_vars[1]
  independent_vars <- formula_vars[-1]

  if (!all(c("id", "time", dependent_var, independent_vars) %in% colnames(data))) {
    stop("Data must contain 'id', 'time', dependent, and independent variables.")
  }

  Y <- reshape_to_matrix(data, "id", "time", dependent_var)
  X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
  X_combined <- do.call(rbind, X_list)
  Z <- rbind(Y, X_combined)

  N <- nrow(Y)
  T <- ncol(Y)
  s <- svd(Z)
  R <- min(N, T)
  Rhat <- which.max(s$d[1:R - 1] / s$d[2:R])
  Mv <- diag(T) - s$v[, Rhat] %*% t(s$v[, Rhat])
  tY <- c(Y %*% Mv)
  tX_combined <- do.call(cbind, lapply(X_list, function(X) c(X %*% Mv)))
  new_data <- data.frame(tY, tX_combined, data[, "id"], data[, "time"])
  colnames(new_data) <- c(formula_vars, "id", "time")

  list(res = plm(formula, data = new_data, index = c("id", "time"), model = "pooling"))
}

# Two-Way BLM Estimator
est_BLM_two_way <- function(formula, data) {
  formula_vars <- all.vars(formula)
  dependent_var <- formula_vars[1]
  independent_vars <- formula_vars[-1]

  if (!all(c("id", "time", dependent_var, independent_vars) %in% colnames(data))) {
    stop("Data must contain 'id', 'time', dependent, and independent variables.")
  }

  Y <- reshape_to_matrix(data, "id", "time", dependent_var)
  X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
  X_combined <- do.call(cbind, X_list)

  N <- nrow(Y)
  T <- ncol(Y)

  clusteri <- cluster_general(Y, X_list, N, T, type = "long")
  G <- clusteri$clusters
  klong <- clusteri$res

  clustert <- cluster_general(Y, X_list, N, T, type = "tall")
  C <- clustert$clusters
  ktall <- clustert$res

  D <- matrix(0, N * T, G * C)

  for (i in seq_len(G)) {
    for (j in seq_len(C)) {
      Du <- as.numeric(klong$cluster == i)
      Dv <- as.numeric(ktall$cluster == j)
      D[, i + (j - 1) * G] <- c(Du %*% t(Dv))
    }
  }
  M <- diag(N * T) - D %*% solve(t(D) %*% D) %*% t(D)

  tY <- M %*% c(Y)
  tX_combined <- do.call(cbind, lapply(X_list, function(X) M %*% c(X)))

  new_data <- data.frame(tY, tX_combined, data[, "id"], data[, "time"])
  colnames(new_data) <- c(formula_vars, "id", "time")

  list(res = plm(formula, data = new_data, index = c("id", "time"), model = "pooling"), df = G * C)
}

# One-Way BLM Estimator
est_BLM_one_way <- function(formula, data) {
  formula_vars <- all.vars(formula)
  dependent_var <- formula_vars[1]
  independent_vars <- formula_vars[-1]

  if (!all(c("id", "time", dependent_var, independent_vars) %in% colnames(data))) {
    stop("Data must contain 'id', 'time', dependent, and independent variables.")
  }

  Y <- reshape_to_matrix(data, "id", "time", dependent_var)
  X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
  X_combined <- do.call(cbind, X_list)

  N <- nrow(Y)
  T <- ncol(Y)

  clusteri <- cluster_general(Y, X_list, N, T, type = "long")
  G <- clusteri$clusters
  klong <- clusteri$res

  D <- matrix(0, N * T, G)

  for (i in seq_len(G)) {
    Du <- as.numeric(klong$cluster == i)
    Dv <- rep(1, T)
    D[, i] <- c(Du %*% t(Dv))
  }

  M <- diag(N * T) - D %*% solve(t(D) %*% D) %*% t(D)

  tY <- M %*% c(Y)
  tX_combined <- do.call(cbind, lapply(X_list, function(X) M %*% c(X)))

  new_data <- data.frame(tY, tX_combined, data[, "id"], data[, "time"])
  colnames(new_data) <- c(formula_vars, "id", "time")

  list(res = plm(formula, data = new_data, index = c("id", "time"), model = "pooling"), df = G)
}

