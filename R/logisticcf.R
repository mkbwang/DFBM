

cf <- function(S, Z, K, lambda, max_iter=1000, tol=1e-4, cap=10, init_A=NULL, init_B=NULL){

  # S: user-item matrix
  # K: number of latent factors
  # lambda: regularization parameter
  # max_iter: maximum number of iterations
  # tol: tolerance
  # cap: limit the maximum of the absolute value in S

  # initialize parameters
  express <- FALSE # express calculation if all entries are observed
  if (mean(Z) == 1) express <- TRUE
  n_obs <- sum(Z)
  full_mask <- Z>0
  N <- nrow(S) # number of rows
  P <- ncol(S) # number of columns
  if (is.null(init_A)){
    A <- matrix(rnorm(N*K, mean=0, sd=0.1), K, N)
    B <- matrix(rnorm(P*K, mean=0, sd=0.1), K, P)
    B[1, ] <- colSums(S*Z)/colSums(Z)
  } else{
    A <- init_A
    B <- init_B
  }
  iter <- 0
  S_lowrank <- t(A) %*% B
  mse <- sum((S[full_mask] - S_lowrank[full_mask])^2) / n_obs
  mse_trace <- rep(0, max_iter)
  # main loop
  while(iter < max_iter){
    if (express){
      new_A <- Matrix::solve(B%*%t(B) + lambda*diag(K), B %*% t(S))
      new_B <- Matrix::solve(new_A%*% t(new_A) + lambda*diag(K), new_A %*% S)
    } else{
      new_A <- A
      new_B <- B
      for(i in 1:N){
        mask <- Z[i, ] > 0
        B_filter <- B[, mask]
        S_filter <- S[i, mask]
        left_mat <- B_filter %*% t(B_filter) + lambda*diag(K)
        right_mat <- B_filter %*% S_filter
        new_A[, i] <- Matrix::solve(left_mat, right_mat)
      }
      for (j in 1:P){
        mask <- Z[, j] > 0
        A_filter <- new_A[, mask]
        S_filter <- S[mask, j]
        left_mat <- A_filter %*% t(A_filter) + lambda * diag(K)
        right_mat <- A_filter %*% S_filter
        new_B[, j] <- Matrix::solve(left_mat, right_mat)
      }
    }
    iter <- iter + 1
    new_S_lowrank <- t(new_A) %*% new_B
    A <- new_A
    B <- new_B
    new_mse <- sum((S[full_mask] - new_S_lowrank[full_mask])^2) / n_obs
    mse_trace[iter] <- new_mse
    if(mse - new_mse < tol){
      mse <- new_mse
      S_lowrank <- new_S_lowrank
      break
    } else{
      mse <- new_mse
      S_lowrank <- new_S_lowrank
    }
  }
  S_lowrank[S_lowrank > cap] <- cap
  S_lowrank[S_lowrank < -cap] <- -cap

  return(list(A=A, B=B, S_lowrank=S_lowrank, n_iter = iter, mse_trace=mse_trace))
}

logisticcf <- function(X, Z, K, lambda, max_iter=1000, tol=1e-5, cap=10){
  # X: user-item matrix
  # K: number of latent factors
  # lambda: regularization parameter
  # max_iter: maximum number of iterations
  # tol: tolerance
  # cap: limit the maximum of the absolute value in S

  # initialize parameters
  N <- nrow(X) # number of rows
  P <- ncol(X) # number of columns
  full_mask <- Z > 0
  n_obs <- sum(Z)
  W <- 2*X - 1

  init_values <- cf(S=4*W, Z=Z, K=K, lambda=lambda, max_iter=max_iter, tol=tol, cap=cap)
  A <- init_values$A
  B <- init_values$B
  S <- init_values$S_lowrank
  pi <- 1/(1 + exp(-S))
  llk <- sum(X[full_mask]*log(pi[full_mask]) + (1-X[full_mask])*log(1-pi[full_mask])) / n_obs
  loss <- rep(0, max_iter)
  iter <- 0
  while(iter < max_iter){
    S_star <- S + 4*W/(1+exp(W*S))
    new_values <- cf(S=S_star, Z=Z, K=K, lambda=lambda, max_iter=max_iter, tol=tol, cap=cap,
                     init_A=A, init_B=B)
    A <- new_values$A
    B <- new_values$B
    S <- new_values$S_lowrank
    pi <- 1/(1 + exp(-S))
    new_llk <- sum(X[full_mask]*log(pi[full_mask]) + (1-X[full_mask])*log(1-pi[full_mask])) / n_obs
    iter <- iter + 1
    loss[iter] <- -new_llk
    if (new_llk - llk < tol){
      llk <- new_llk
      break
    } else{
      llk <- new_llk
    }
  }
  return(list(A=A, B=B, pi=pi, loss_trace=loss[1:iter]))
}
