source("experiment/visualization_utils.R")

count_mat <- read.table("experiment/tabula_muris/ZINBwave/simulated_counts.tsv",
                        header=TRUE, row.names=1, sep="\t")


library(Matrix)

slice_0 <- 1*(count_mat > 0)
mask <- matrix(1, nrow=nrow(count_mat),
                       ncol=ncol(count_mat))
nonzero_indices <- which(mask > 0, arr.ind=T)


train_index <- sample(nrow(nonzero_indices), nrow(nonzero_indices)*0.8)
train_indices <- nonzero_indices[train_index, ]
validation_indices <- nonzero_indices[-train_index, ]
mask_train <- mask
mask_validation <- mask

mask_train[validation_indices] <- 0
mask_validation[train_indices] <- 0

slice_0_train <- slice_0
slice_0_train[validation_indices] <- NA


# fit with cxx implementation of logisticcf
begin <- proc.time()
result_0_cxx <- logisticcf(slice_0, mask_train, 10, 0.1, 1000, 1e-5)
cxx_validation_probs <- result_0_cxx$pi[validation_indices]
end <- proc.time()
end-begin

# fit with R implementation of logisticcf
begin <- proc.time()
result_0_r <- logisticcfR(slice_0, mask_train, K=1, lambda=0.1, tol=2e-5)
r_validation_probs <- result_0_r$pi[validation_indices]
end <- proc.time()
end-begin

truth_validation <- slice_0[validation_indices]

library(pROC)

roc(truth_validation, r_validation_probs)
mean(truth_validation*log(r_validation_probs)+
       (1-truth_validation)*log(1-r_validation_probs))

roc(truth_validation, cxx_validation_probs)
mean(truth_validation*log(cxx_validation_probs)+
       (1-truth_validation)*log(1-cxx_validation_probs))

library(logisticPCA)

# fit with logisticSVD
begin <- proc.time()
lsvd_fit <- logisticSVD(x=slice_0_train, k=2, conv_criteria=1e-5)
lsvd_probs <- fitted(lsvd_fit, type="response")
end <- proc.time()
end - begin

lsvd_validation_probs <- lsvd_probs[validation_indices]
roc(truth_validation, lsvd_validation_probs)
mean(truth_validation*log(lsvd_validation_probs)+
       (1-truth_validation)*log(1-lsvd_validation_probs))



# test an extremely sparse case

mask <- matrix(0, nrow=10, ncol=10)
indices <- cbind(sample(10), sample(10))
mask[indices] <- 1

values <- matrix(rbinom(n=100, size=1, prob=0.6),
                  nrow=10, ncol=10)
observed_values <- values * mask

result_r <- logisticcfR(observed_values, mask, K=3, lambda=0.01, tol=1e-5)
result_cxx <- logisticcf(observed_values, mask, 1, 0.1, 1000, 1e-5)



