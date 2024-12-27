source("experiment/visualization_utils.R")

count_mat <- read.table("experiment/tabula_muris/ZINBwave/simulated_counts.tsv",
                        header=TRUE, row.names=1, sep="\t")


library(Matrix)

slice_0 <- 1*(count_mat > 0)
nonzero_indices <- which(slice_0 > 0, arr.ind=T)


train_index <- sample(nrow(nonzero_indices), nrow(nonzero_indices)*0.8)
train_indices <- nonzero_indices[train_index, ]
validation_indices <- nonzero_indices[-train_index, ]
slice_0_train <- slice_0
slice_0_validation <- slice_0
slice_0_train[validation_indices] <- 0
slice_0_validation[train_indices] <- 0

slice_5 <- 1*(count_mat > 5)
slice_10 <- 1*(count_mat > 10)

initial_mask <- matrix(1, nrow=nrow(count_mat),
                       ncol=ncol(count_mat))


library(logisticPCA)
# begin <- proc.time()
# lpca_fit <- logisticPCA(x=slice_0, k=4, m=5)
# end <- proc.time()
# end - begin
# lpca_probs <- fitted(lpca_fit, type="response")

slice_10_0_train <- matrix(NA, nrow=nrow(slice_10), ncol=ncol(slice_10))
slice_10_0_train[train_indices] <- slice_10[train_indices]

slice_5_0_train <- matrix(NA, nrow=nrow(slice_5), ncol=ncol(slice_5))
slice_5_0_train[train_indices] <- slice_5[train_indices]


begin <- proc.time()
lsvd_fit <- logisticSVD(x=slice_5_0_train, k=3)
lsvd_probs <- fitted(lsvd_fit, type="response")
end <- proc.time()
end - begin
lsvd_probs_validation <- lsvd_probs[validation_indices]


begin <- proc.time()
lcf_fit <- logisticcf(X=slice_5, Z=slice_0_train, K=3, lambda=0.09, tol=1e-5)
lcf_probs <- lcf_fit$pi
lcf_probs_validation <- lcf_probs[validation_indices]
end <- proc.time()
end - begin

truth_validation <- slice_5[validation_indices]

library(pROC)

roc(truth_validation, lsvd_probs_validation)
lsvd_validation_mean_llk <- mean(truth_validation*log(lsvd_probs_validation)+
                                   (1-truth_validation)*log(1-lsvd_probs_validation))

roc(truth_validation, lcf_probs_validation)
lcf_validation_mean_llk <- mean(truth_validation*log(lcf_probs_validation)+
                                  (1-truth_validation)*log(1-lcf_probs_validation))

#
# begin <- proc.time()
# nbmf_fit <- nbmf(Y= t(slice_0), Z=t(initial_mask), k=3, tol=1e-5)
# nbmf_probs = t(nbmf_fit$pi)
# end <- proc.time()
# end - begin
