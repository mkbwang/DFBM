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
slice_30 <- 1*(count_mat > 30)
slice_100 <- 1*(count_mat > 100)
slice_200 <- 1*(count_mat > 200)


initial_mask <- matrix(1, nrow=nrow(count_mat),
                       ncol=ncol(count_mat))

# fit with cxx implementation of logisticcf
begin <- proc.time()
result_cxx <- logisticcf(slice_5, slice_0, 4, 10, 1000, 1e-5, 10)
# cxx_validation_probs <- result_0_cxx$pi[validation_indices]
end <- proc.time()
end-begin

# fit with R implementation of logisticcf
begin <- proc.time()
result_0_r <- logisticcfR(slice_5, slice_0, K=4, lambda=10, tol=1e-5)
# r_validation_probs <- result_0_r$pi[validation_indices]
end <- proc.time()
end-begin

truth_validation <- slice_5[validation_indices]

library(pROC)

roc(truth_validation, cxx_validation_probs)
roc(truth_validation, r_validation_probs)
mean(truth_validation*log(cxx_validation_probs)+
                                   (1-truth_validation)*log(1-cxx_validation_probs))
mean(truth_validation*log(r_validation_probs)+
       (1-truth_validation)*log(1-r_validation_probs))


library(logisticPCA)

# fit with logisticSVD
begin <- proc.time()
lsvd_fit <- logisticSVD(x=slice_0, k=4, conv_criteria=1e-5)
lsvd_probs <- fitted(lsvd_fit, type="response")
end <- proc.time()
end - begin


# train and provide validation loss
begin <- proc.time()
result_lcf <- train_lcf(train_mask = slice_0_train,
                        validation_mask = slice_0_validation,
                        observation = slice_5,
                        K=4, lambda=1, max_iter=1000, tol=1e-5, cap=10)
end <- proc.time()
end - begin


