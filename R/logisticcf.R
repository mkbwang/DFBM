



log1exp <- function(x){
  output <- x
  output[x<20] <- log(1+exp(x[x<20]))
  return(output)
}


penalized_loss <- function(X, Z, A, B, lambda){
  # X: binary matrix, N*P
  # Z: binary mask, N*P
  # A: latent matrix A, K*N
  # B: latent matrix B, K*P
  # lambda: regularization parameter
  S <- t(A) %*% B
  logloss <- sum(-Z*(X*S - log1exp(S)))
  ridge_penalty <- lambda * (sum(A[-1, ]^2) + sum(B[-1, ]^2))
  return((logloss + ridge_penalty)/sum(Z))
}


#' @useDynLib NDBEC
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
logisticcfR <- function(X, Z, K, lambda, max_iter=1000, tol=1e-5){
  # X: user-item matrix
  # K: number of latent factors
  # lambda: regularization parameter
  # max_iter: maximum number of iterations
  # tol: tolerance

  # initialize parameters
  N <- nrow(X) # number of rows
  P <- ncol(X) # number of columns
  W <- 2*X - 1
  XZ <-  X * Z

  # initialize A
  A <- matrix(rnorm(N*K, mean=0, sd=0.1), K, N)
  A[1, ] <- rowSums(W*Z)/rowSums(Z)
  B <- matrix(rnorm(P*K, mean=0, sd=0.1), K, P)
  B[1, ] <- colSums(W*Z)/colSums(Z)
  S <- t(A) %*% B

  A_penalty <- matrix(lambda, nrow=K, ncol=N)
  A_penalty[1, ] <- 0
  B_penalty <- matrix(lambda, nrow=K, ncol=P)
  B_penalty[1, ] <- 0

  loss <- penalized_loss(X, Z, A, B, lambda)
  loss_trace <- rep(0, max_iter+1)
  loss_trace[1] <- loss
  iter <- 1
  while(iter <= max_iter){
    # update A
    pi <- 1/(1+exp(-S))
    g_star <- -B %*% (t(XZ)) + B %*% (t(Z)*t(pi)) + 2 * A_penalty * A
    f_star <- 0.25 * (B*B) %*% t(Z) + 2*A_penalty
    A <- A - g_star/(2*f_star)

    # update B
    S_dagger <- t(A) %*% B
    pi_dagger <- 1/(1+exp(-S_dagger))
    g_dagger <- -A %*% XZ + A %*% (Z*pi_dagger) + 2* B_penalty * B
    f_dagger <- 0.25 * (A*A) %*% Z + 2*B_penalty
    B <- B - g_dagger/(2*f_dagger)

    S <- t(A) %*% B
    new_loss <- penalized_loss(X, Z, A, B, lambda)
    iter <- iter  + 1
    loss_trace[iter] <- new_loss
    if (loss - new_loss < tol){
      loss <- new_loss
      break
    } else{
      loss <- new_loss
    }
  }
  pi <- 1/(1 + exp(-S))
  return(list(A=A, B=B, pi=pi, loss_trace=loss_trace[1:iter]))
}
