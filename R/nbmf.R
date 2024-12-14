




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

  loglik <- sum(Z*(Y*log(pi) + (1-Y)*log(1-pi)))
  # print(sprintf("Maximum value in A: %f", max(A)))
  # print(sprintf("Minimum value in A: %f", min(A)))
  # print(sprintf("Maximum value in B: %f", max(B)))
  # print(sprintf("Minimum value in B: %f", min(B)))
  # print(sprintf("Maximum value in pi: %f", max(pi)))
  # print(sprintf("Minimum value in pi: %f", min(pi)))
  # print(sprintf("log-likelihood: %f", loglik))
  # browser()
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
    # new_pi[new_pi > 1-eps] <- 1-eps
    # new_pi[new_pi < eps] <- eps
    new_loglik <- sum(Z*(Y*log(new_pi) + (1-Y)*log(1-new_pi)))

    # print(sprintf("Maximum value in A: %f", max(new_A)))
    # print(sprintf("Minimum value in A: %f", min(new_A)))
    # print(sprintf("Maximum change in A: %f", max(abs(new_A - A))))
    # print(sprintf("Maximum value in B: %f", max(new_B)))
    # print(sprintf("Minimum value in B: %f", min(new_B)))
    # print(sprintf("Maximum change in B: %f", max(abs(new_B - B))))
    # print(sprintf("Maximum value in pi: %f", max(new_pi)))
    # print(sprintf("Minimum value in pi: %f", min(new_pi)))
    # print(sprintf("Maximum change in pi: %f", max(abs(new_pi - pi))))
    # print(sprintf("log-likelihood: %f", new_loglik))
    # print(sprintf("Change in log-likelihood: %f", new_loglik - loglik))
    # browser()
    A <- new_A
    B <- new_B
    pi <- new_pi
    # check convergence
    if (abs(new_loglik - loglik) / sum(Z)  < tol){
      break
    } else{
      loglik <- new_loglik
    }

  }
  pi_all[rows_include, columns_include] <- pi
  converge <- ifelse(iter < max_iter, TRUE, FALSE)
  return(list(rows = rows_include, columns=columns_include,
              A=A, B=B, pi=pi_all, loglik = loglik, iter=iter, converge=converge))
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
