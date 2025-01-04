

#' Train_LCF
#' @description
#' training logistic collaborative filtering and test it on a validation set
#'
#' @useDynLib NDBEC, .registration=TRUE
#' @importFrom pROC auc
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @param train_mask a binary matrix indicating the training set
#' @param validation_mask a binary matrix indicating the validation set
#' @param observation a matrix of observed values, 0 or 1
#' @param K dimension of latent space
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations for collaborative filtering
#' @param tol tolerance for collaborative filtering, change of logit loss between iterations
#' @param cap limit the maximum of the absolute value of estimated logit link
#' @returns logistic collaborative filtering results and performance on validation set
#' @export
train_lcf <- function(train_mask, validation_mask, observation,
                      K,  lambda, max_iter=1000, tol=1e-5, cap=10){

  fitted_result <- logisticcf(observation, train_mask, K, lambda, max_iter, tol, cap)
  A_mat <- fitted_result$A
  B_mat <- fitted_result$B
  pi_mat <- fitted_result$pi
  loss_trace <- fitted_result$loss_trace
  validation_indices <- which(validation_mask > 0, arr.ind=T)
  pi_validation <- pi_mat[validation_indices]
  truth_validation <- observation[validation_indices]

  auc <- pROC::auc(truth_validation, pi_validation)
  llk_validation <- mean(truth_validation*log(pi_validation)+
                                    (1-truth_validation)*log(1-pi_validation))

  return(list(A=A_mat, B=B_mat, pi=pi_mat, loss_trace=loss_trace,
              auc_validation=auc, llk_validation=llk_validation))

}


cfR <- function(S, Z, K, lambda, max_iter=1000, tol=1e-4, cap=10, init_A=NULL, init_B=NULL){

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
