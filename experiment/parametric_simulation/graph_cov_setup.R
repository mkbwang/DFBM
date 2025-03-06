

rm(list=ls())


# first generate a graph (erdos renyi)
dimension <- 10
mygraph <- matrix(0, nrow=dimension, ncol=dimension)

upper_indices <- which(upper.tri(mygraph, diag=F))
nz_indices <- sample(upper_indices, size=0.2*length(upper_indices))
mygraph[nz_indices] <- 1
mygraph <- mygraph + t(mygraph)

# then generate off diagonal points of the precision matrix
prec_mat <- matrix(0, nrow=dimension, ncol=dimension)
lbound <- 2
ubound <- 3
pos_indices <- sample(nz_indices, size=0.5*length(nz_indices))
neg_indices <- setdiff(nz_indices, pos_indices)
prec_mat[pos_indices] <- runif(length(pos_indices), min=lbound, max=ubound)
prec_mat[neg_indices] <- runif(length(neg_indices), min=-ubound, max=-lbound)
prec_mat <- prec_mat + t(prec_mat)


# the choice of diagonal values are crucial
# since they guarantees that the matrix is positive definite and tunes the condition number
## initial values
diag(prec_mat) <- max(abs(prec_mat)) * 2


library(RSpectra)
accurate_kappa <- function(symmat){
  largest <- eigs_sym(symmat, k=1, which="LM")$values
  smallest <- eigs_sym(symmat, k=1, sigma=0)$values
  largest/smallest
}

# function that searches for
binSearchCond <- function(Theta, condTheta, numBinSearch=100, epsBin=0.01) {
  # Internal function that determines the constant in the diagonal to satisfy the
  # condition constraint on the Precision/Covariance matrix

  n <- nrow(Theta)
  currCondTheta <- accurate_kappa(Theta)
  if (currCondTheta < condTheta) {
    # Max entry in the diagonal (lower bound)
    currLB   <- -max(diag(Theta))
    stepSize <- currLB+.Machine$double.eps

    while (currCondTheta < condTheta) {
      currCondTheta <- accurate_kappa(Theta+stepSize*diag(n))
      stepSize      <- stepSize/2
    }
    currUB <- stepSize
  } else {
    currLB <- 0
    stepSize = 0.1

    while (currCondTheta > condTheta) {
      currCondTheta <- accurate_kappa(Theta + stepSize*diag(n))
      stepSize      <- 2*stepSize
    }
    currUB <- stepSize
  }

  for (i in 1:numBinSearch) {
    diagConst <- (currUB+currLB)/2
    currCondTheta <- accurate_kappa(Theta+diagConst*diag(n))

    if (currCondTheta < condTheta) currUB <- diagConst
    else currLB <- diagConst

    if (abs(currCondTheta-condTheta)<epsBin) break
  }
  diagConst
}

adjustment <- binSearchCond(Theta=prec_mat, condTheta=10)
prec_mat <- prec_mat + adjustment*diag(dimension)
## check condition number
accurate_kappa(prec_mat)


# take inverse to get covariance matrix
covar_mat <- chol2inv(chol(prec_mat))
corr_mat <- cov2cor(covar_mat)

