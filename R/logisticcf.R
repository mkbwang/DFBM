

cf <- function(S, K, lambda, max_iter=1000, tol=1e-4){

  # S: user-item matrix
  # K: number of latent factors
  # lambda: regularization parameter
  # max_iter: maximum number of iterations
  # tol: tolerance

  # initialize parameters
  N <- nrow(S) # number of rows
  P <- ncol(S) # number of columns
  A <- matrix(rnorm(N*K, mean=0, sd=0.1), K, N)
  B <- matrix(rnorm(P*K, mean=0, sd=0.1), K, P)
  B[1, ] <- colMeans(S)
  iter <- 0
  S_lowrank <- t(A) %*% B
  mse <- sum((S - S_lowrank)^2)/(N*P)
  mse_trace <- rep(0, max_iter)
  # main loop
  while(iter < max_iter){
    # first update A
    new_A = Matrix::solve(B%*%t(B) + lambda*diag(K), B %*% t(S))
    new_B = Matrix::solve(new_A%*% t(new_A) + lambda*diag(K), new_A %*% S)
    iter <- iter + 1
    new_S_lowrank <- t(new_A) %*% new_B
    A <- new_A
    B <- new_B
    new_mse <- sum((S - new_S_lowrank)^2)/(N*P)
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
  return(list(A=A, B=B, S_lowrank=S_lowrank, n_iter = iter))
}

logisticcf <- function(X, K, lambda, max_iter=1000, tol=1e-4){
  # X: user-item matrix
  # K: number of latent factors
  # lambda: regularization parameter
  # max_iter: maximum number of iterations
  # tol: tolerance

  # initialize parameters
  N <- nrow(X) # number of rows
  P <- ncol(X) # number of columns
  W <- 2*X - 1

  init_values <- cf(S=4*W, K=K, lambda=lambda, max_iter=max_iter, tol=tol)
  A <- init_values$A
  B <- init_values$B
  S <- t(A) %*% B
  pi <- 1/(1 + exp(-S))
  llk <- mean(X*log(pi) + (1-X)*log(1-pi))
  loss <- rep(0, max_iter)
  iter <- 0
  while(iter < max_iter){
    S_star <- S + 4*W/(1+exp(W*S))
    new_values <- cf(S=S_star, K=K, lambda=lambda, max_iter=max_iter, tol=tol)
    A <- new_values$A
    B <- new_values$B
    S <- t(A) %*% B
    pi <- 1/(1 + exp(-S))
    new_llk <- mean(X*log(pi) + (1-X)*log(1-pi))
    iter <- iter + 1
    loss[iter] <- -new_llk
    if (new_llk - llk < tol){
      llk <- new_llk
      break
    } else{
      llk <- new_llk
    }
  }
  return(list(A=A, B=B, pi = pi, loss_trace=loss[1:iter]))
}
