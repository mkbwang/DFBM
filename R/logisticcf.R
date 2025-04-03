
#' log one plus exponential
#'
#' @param x a vector
#' @returns log(1+exp(x))
#' @details
#' This function is used to avoid numerical issues when x is too large
#' @export
log1exp <- function(x){
  output <- x
  output[x<20] <- log(1+exp(x[x<20]))
  return(output)
}


#' Cross entropy loss plus ridge penalty
#'
#' @param X observed binary matrix
#' @param Z binary mask indicating whether an entry is included in the loss
#' @param A factor matrix A
#' @param B factor matrix B expit(A^T*B) is the denoised probability matrix
#' @param lambda ridge penalty parameter
#' @return value of loss
#' @export
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



#' Logistic collaborative filtering
#'
#' @param X observed binary matrix
#' @param Z binary mask indicating whether an entry should be considered when fitting logisticcf
#' @param K the maximum number of ranks of factorized matrices
#' @param lambda ridge penalty parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence limit of the loss
#' @importFrom stats rnorm
#' @return a list containting the factorized functions and the  denoised probabilities
#' @export
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
    A[f_star > 0] <- A[f_star > 0] - g_star[f_star > 0]/(2*f_star[f_star > 0])
    # update B
    S_dagger <- t(A) %*% B
    pi_dagger <- 1/(1+exp(-S_dagger))
    g_dagger <- -A %*% XZ + A %*% (Z*pi_dagger) + 2* B_penalty * B
    f_dagger <- 0.25 * (A*A) %*% Z + 2*B_penalty
    B[f_dagger > 0] <- B[f_dagger > 0] - g_dagger[f_dagger > 0]/(2*f_dagger[f_dagger > 0])
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


#' Cross validation to select optimal rank of logisticcf
#' @param X observed binary matrix
#' @param Z binary mask indicating whether an entry should be considered when fitting logisticcf
#' @param max_K maximum number of ranks to try for matrix factorization
#' @param lambdas ridge penalty parameters
#' @param ncores number of cores for parallel computing, default 1
#' @returns optimal choice of factorizedmatrix A and B, denoised probabilities, selected rank
#'
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom parallelly availableCores
#' @importFrom foreach foreach %dopar%
#' @export
cv.logisticcfR <- function(X, Z, max_K=10, lambdas=c(0.01, 0.1, 1), ncores=1){

  dim1 <- nrow(X)
  dim2 <- ncol(X)

  # set up training and validation sets
  Z_train <- matrix(0, nrow=dim1, ncol=dim2)
  Z_validation <- matrix(0, nrow=dim1, ncol=dim2)

  # split the entries into training and validation
  nonzero_mask_indices <- which(Z > 0)
  train_indices <- sample(nonzero_mask_indices,
                          size=0.8*length(nonzero_mask_indices))
  validation_indices <- setdiff(nonzero_mask_indices, train_indices)
  Z_train[train_indices] <- 1
  Z_validation[validation_indices] <- 1

  # filter out rows and columns that have completely missing entries
  row_sums <- rowSums(Z_train)
  col_sums <- colSums(Z_train)
  X_subset <- X
  if (any(row_sums== 0) | any(col_sums == 0)){
    # some rows and columns need to be filtered
    selected_rows <- which(row_sums > 0)
    selected_cols <- which(col_sums > 0)
    X_subset <- X[selected_rows, selected_cols, drop=FALSE]
    Z_train <- Z_train[selected_rows, selected_cols, drop=FALSE]
    Z_validation <- Z_validation[selected_rows, selected_cols, drop=FALSE]
  }

  upper_bound_K <- ceiling(sum(Z_train)/(nrow(Z_train) + ncol(Z_train)))
  max_K <- min(max_K, upper_bound_K)

  validation_losses <- matrix(0, nrow=max_K, ncol=length(lambdas))

  rank_1_fit <- logisticcfR(X_subset, Z_train, K=1, lambda=0)
  validation_losses[1, ] <- penalized_loss(X_subset, Z_validation,
                                           rank_1_fit$A, rank_1_fit$B,
                                           lambda=0)
  # (sprintf("Rank 1 validation loss is: %f", validation_losses[1, 1]))
  selected_K <- 1
  selected_lambda <- 0
  numCores <- min(availableCores(), ncores)
  if (max_K > 1){
    for (k in 2:max_K){
      # print(sprintf("Fitting rank %d models", k))
      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      clusterExport(cl, varlist=c("logisticcfR", "penalized_loss", "log1exp"))
      j <- 1
      loss_vec <- foreach(j=1:length(lambdas),
                          .combine=c) %dopar%{
                            fit <- logisticcfR(X_subset, Z_train, K=k, lambda=lambdas[j])
                            validation_losses[k, j] <- penalized_loss(X_subset, Z_validation,
                                                                      fit$A, fit$B, lambda=0)
                          }
      validation_losses[k, ] <- loss_vec
      stopCluster(cl)
      # print(sprintf("Minimum Rank %d validation loss is: %f", k, min(validation_losses[k, ])))
      if (min(validation_losses[k, ]) < min(validation_losses[k-1, ])){
        selected_K <- k
        selected_lambda <- lambdas[which.min(validation_losses[k, ])]
      } else{ # end early if the validation loss begins to increase
        break
      }
    }
  }
  # print("Fitting the final model...")
  # refit a final model with the selected K and lambda
  final_model <- logisticcfR(X, Z, K=selected_K, lambda=selected_lambda)


  return(list(A=final_model$A, B=final_model$B, pi=final_model$pi,
              selected_K = selected_K, selected_lambda=selected_lambda,
              validation_losses=validation_losses[1:selected_K, ]))

}
