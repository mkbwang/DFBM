

negloglik <- function(Y, Z, pi, avg = F){

  sum_loglik <- sum(Z*(Y*log(pi) + (1-Y)*log(1-pi)))
  if (avg){
    avg_loglik <- sum_loglik / sum(Z)
  }
  if (avg){
    return(-avg_loglik)
  } else{
    return(-sum_loglik)
  }

}


nbmf <- function(Y, Z, k=2, max_iter=500, eps=1e-6, tol=1e-5){

  # Y: binary matrices indicating whether count matrix X is larger than a threshold
  # Z: binary matrices indicating missingness (0). A missing entry is not included in the factorization
  # k: rank of the factorization
  # max_iter: maximum number of iterations
  # tol: tolerance for convergence
  # returns: factorized matrices A and B

  pi_all <-matrix(0, nrow=nrow(Y), ncol=ncol(Y))

  # first get rid of columns and rows where Z is completely zero
  column_sum <- colSums(Z)
  columns_include <- which(column_sum > 0)
  row_sum <- rowSums(Z)
  rows_include <- which(row_sum > 0)
  Y <- Y[rows_include, columns_include]
  Z <- Z[rows_include, columns_include]
  ZY <- Z*Y

  N <- nrow(Y)
  P <- ncol(Y)
  Lambda <- matrix(1, nrow=k, ncol=N) %*% Z

  ## check if the rank is feasible given the available values for estimation
  num_params <- N*k + k*P - P
  num_obs <- sum(Z) # number of observed values
  if (num_params > num_obs){
    stop("The rank is too high for parameters to be estimated")
  }

  # initialize A and B
  observed_prob <- sum(ZY)/sum(Z)
  A <- matrix(rbeta(N*k, shape1=observed_prob, shape2=1-observed_prob),
              nrow=N, ncol=k)
  A[A > 1-eps] <- 1-eps
  A[A < eps] <- eps

  B <- matrix(runif(k*P), nrow=k, ncol=P)
  B <- B / (matrix(1, nrow=k, ncol=k) %*% B)

  pi <- A %*% B

  loss <- negloglik(Y, Z, pi)
  # iterate until convergence
  iter <- 0
  while (iter < 500){

    iter <- iter + 1
    # update A
    C <- A * ((ZY/pi) %*% t(B))
    H <- (1-A) * (((Z-ZY)/(1-pi)) %*% t(B))
    new_A <- C/(C+H)

    new_A[new_A > 1-eps] <- 1-eps
    new_A[new_A < eps] <- eps

    pi_dagger <- new_A %*% B
    # update B
    new_B <- B * (t(new_A) %*% (ZY/pi_dagger)) / Lambda +
      B * ((1-t(new_A)) %*% ((Z-ZY)/(1-pi_dagger))) / Lambda

    new_B[new_B > 1-eps] <- 1-eps
    new_B[new_B < eps] <- eps
    new_B <- new_B / (matrix(1, nrow=k, ncol=k) %*% new_B)

    new_pi <- new_A %*% new_B
    new_loss <- negloglik(Y, Z, new_pi)

    A <- new_A
    B <- new_B
    pi <- new_pi
    # check convergence
    if (abs(new_loss - loss) / sum(Z)  < tol){
      break
    } else{
      loss <- new_loss
    }
  }
  AIC <- 2 * loss + 2 * num_params
  pi_all[rows_include, columns_include] <- pi
  converge <- ifelse(iter < max_iter, TRUE, FALSE)
  return(list(rows = rows_include, columns=columns_include,
              A=A, B=B, pi=pi_all, loss = loss, AIC=AIC,
              iter=iter, converge=converge))
}



cv.nbmf <- function(Y, Z, max_k=10, fold=5, max_iter=500,
                    eps=1e-6, tol=1e-5, seed = 1){

  set.seed(seed)
  mask_coordinates <- which(Z > 0, arr.ind=T)
  nrows <- length(unique(mask_coordinates[, 1]))
  ncols <- length(unique(mask_coordinates[, 2]))

  # confirm the largest rank possible
  num_obs <- sum(Z)
  maximum_feasible_rank <- floor( (num_obs+ ncols)/(nrows+ncols) )
  if (max_k > maximum_feasible_rank){
    max_k <- maximum_feasible_rank
  }

  # create the CV folds
  indices <- sample(1:num_obs)
  cvblocks <- list()
  begin_index <- as.integer(num_obs * seq(0, fold-1) / fold)+1
  for (j in 1:fold){
    if (j == fold){
      cvblocks[[j]] <- indices[begin_index[j]:num_obs]
    } else{
      cvblocks[[j]] <- indices[begin_index[j]:(begin_index[j+1]-1)]
    }
  }

  avg_loss <- rep(0, max_k)
  for (rank in 1:max_k){

    print(sprintf("Rank: %d", rank))
    test_losses <- rep(0, fold)

    for (j in 1:fold){

      print(sprintf("Fold number: %d", j))
      train_mask <- sparseMatrix(i = mask_coordinates[unlist(cvblocks[-j]), 1],
                                 j = mask_coordinates[unlist(cvblocks[-j]), 2],
                                 x = rep(1, length(unlist(cvblocks[-j]))),
                                 dims = c(dim(Z)),
                                 repr="T")

      test_mask <- sparseMatrix(i = mask_coordinates[unlist(cvblocks[j]), 1],
                                j = mask_coordinates[unlist(cvblocks[j]), 2],
                                x = rep(1, length(unlist(cvblocks[j]))),
                                dims = c(dim(Z)),
                                repr="T")

      selected_rows <- which(Matrix::rowSums(train_mask) > 0)
      selected_columns <- which(Matrix::colSums(train_mask) > 0)

      observation <- Y[selected_rows, selected_columns]
      train_mask <- train_mask[selected_rows, selected_columns]
      test_mask <- test_mask[selected_rows, selected_columns]

      fit_result <- nbmf(Y=as.matrix(observation), Z=as.matrix(train_mask), k=rank,
                         max_iter=max_iter, eps=eps, tol=tol)

      test_loss <- negloglik(Y = as.matrix(observation),
                             Z = as.matrix(test_mask),
                             pi = fit_result$pi, avg=T)
      test_losses[j] <- test_loss
    }
    avg_loss[rank] <- mean(test_losses)
    print(test_losses)
  }

  best_rank <- which.min(avg_loss)
  return(list(best_rank=best_rank, loss = avg_loss))

}


# slice_factorization <- function(Y, Z, max_k=10, cv=T){
#
#   # Y: binary matrices indicating whether count matrix X is larger than a threshold
#   # Z: binary matrices indicating missingness (0). A missing entry is not included in the factorization
#   # max_k: maximum rank to try for nonnegative binary matrix factorization
#   # cv: perform cross-validation to choose the rank, if not use AIC
#   # returns: denoised and batch-effect corrected binary matrix
#
#
# }
