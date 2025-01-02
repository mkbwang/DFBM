
rm(list=ls())

model_selection_fit <- function(Y, Z, max_k=6, seed=2024){

  cv_result <- cv.nbmf(Y, Z , max_k=max_k, seed=seed)
  AICs <- rep(0, max_k)
  for (j in 1:max_k){
    print(j)
    result <- nbmf(Y, Z, k=j)
    end <- proc.time()
    AICs[j] <- result$AIC
  }

  best_AIC <- which.min(AICs)
  best_rank <- max(cv_result$best_rank, best_AIC)

  selected_result <- nbmf(Y , Z,
                          k=best_rank)

  output <- list(cv_loss = cv_result$loss, AICs=AICs,
                 fit = selected_result)

  return(output)

}

source("experiment/visualization_utils.R")

count_mat <- read.table("experiment/tabula_muris/ZINBwave/simulated_counts.tsv",
                        header=TRUE, row.names=1, sep="\t")


library(Matrix)
slice_0 <- 1*(count_mat > 0)
plot_heatmap(slice_0, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
initial_mask <- matrix(1, nrow=nrow(count_mat),
                       ncol=ncol(count_mat))


slice_2 <- 1*(count_mat > 2)
slice_2_withholes <- slice_2
slice_2_withholes[slice_0 == 0] <- NA

slice_5 <- 1*(count_mat > 5)
slice_5_withholes <- slice_5
slice_5_withholes[slice_2 == 0] <- NA

slice_10 <- 1*(count_mat > 10)
slice_10_withholes <- slice_10
slice_10_withholes[slice_5 == 0] <- NA

# logisticPCA fit
library(logisticPCA)
begin <- proc.time()
lpca_fit <- logisticPCA(x=slice_10_withholes, k=2, m=5)
end <- proc.time()
end - begin
lpca_probs <- fitted(lpca_fit, type="response")


begin <- proc.time()
lsvd_fit_0 <- logisticSVD(x=slice_0, k=2)
lsvd_probs_0 <- fitted(lsvd_fit_0, type="response")
lsvd_fit_2 <- logisticSVD(x=slice_2, k=2)
lsvd_probs_2 <- fitted(lsvd_fit_2, type="response")
end <- proc.time()
end - begin


# NBMF fit
begin <- proc.time()
nbmf_fit <- nbmf(Y= t(slice_10), Z=t(slice_5), k=3)
nbmf_probs = t(nbmf_fit$pi)
end <- proc.time()
end - begin


library(pROC)
auc(as.vector(slice_10[slice_5 > 0]), as.vector(nbmf_probs[slice_5 > 0]))
auc(as.vector(slice_10[slice_5 > 0]), as.vector(lsvd_probs[slice_5 > 0]))

output <- model_selection_fit(Y = slice_0, Z = initial_mask, max_k=14,
                              seed=2024)
begin <- proc.time()
fit_result <- nbmf(Y = slice_0, Z = initial_mask,
     k=8)
end <- proc.time()

par(mar = rep(2, 4))
plot(output$cv_loss, xlab="", ylab="")
plot(output$AICs, xlab="", ylab="")
smoothed_value <- output$fit$pi
plot_heatmap(smoothed_value, legend="estimated", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
slice_1 <- 1*(count_mat > 1)
slice_2 <- 1*(count_mat > 2)
slice_3 <- 1*(count_mat > 3)
slice_4 <- 1*(count_mat > 4)
slice_5 <- 1*(count_mat > 5)
slice_10 <- 1*(count_mat > 10)

output_0_1 <- model_selection_fit(Y = slice_1, Z = slice_0, max_k=5,
                                  seed=2024)

output_1_2 <- model_selection_fit(Y = slice_2, Z = slice_1, max_k=5,
                                  seed=2024)

output_2_3 <- model_selection_fit(Y = slice_3, Z = slice_2, max_k=5,
                                  seed=2024)

output_0_2 <- model_selection_fit(Y = slice_2, Z = slice_0, max_k=5,
                                   seed=2024)

observed_mat <- slice_2
estimated_mat <- output_0_2$fit$pi
observed_mat[slice_0 == 0] <- NA
estimated_mat[slice_0 == 0] <- NA
plot_heatmap(observed_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
plot_heatmap(estimated_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
par(mar = rep(2, 4))
plot(output_0_2$cv_loss, xlab="", ylab="")
plot(output_0_2$AICs, xlab="", ylab="")


output_0_5 <- model_selection_fit(Y = slice_5, Z = slice_0, max_k=5,
                                  seed=2024)

observed_mat <- slice_5
estimated_mat <- output_0_5$fit$pi
observed_mat[slice_0 == 0] <- NA
estimated_mat[slice_0 == 0] <- NA
plot_heatmap(observed_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
plot_heatmap(estimated_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
par(mar = rep(2, 4))
plot(output_0_5$cv_loss, xlab="", ylab="")
plot(output_0_5$AICs, xlab="", ylab="")


output_2_5 <- model_selection_fit(Y = slice_5, Z = slice_2, max_k=5,
                                  seed=2024)


observed_mat <- slice_5
estimated_mat <- output_2_5$fit$pi
observed_mat[slice_2 == 0] <- NA
estimated_mat[slice_2 == 0] <- NA
plot_heatmap(observed_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
plot_heatmap(estimated_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
par(mar = rep(2, 4))
plot(output_2_5$cv_loss, xlab="", ylab="")
plot(output_2_5$AICs, xlab="", ylab="")


output_0_10 <- model_selection_fit(Y = slice_10, Z = slice_0, max_k=5,
                                  seed=2024)

observed_mat <- slice_10
estimated_mat <- output_0_10$fit$pi
observed_mat[slice_0 == 0] <- NA
estimated_mat[slice_0 == 0] <- NA
plot_heatmap(observed_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
plot_heatmap(estimated_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
par(mar = rep(2, 4))
plot(output_0_10$cv_loss, xlab="", ylab="")
plot(output_0_10$AICs, xlab="", ylab="")


output_5_10 <- model_selection_fit(Y = slice_10, Z = slice_5, max_k=5,
                                  seed=2024)

observed_mat <- slice_10
estimated_mat <- output_5_10$fit$pi
observed_mat[slice_5 == 0] <- NA
estimated_mat[slice_5 == 0] <- NA
plot_heatmap(observed_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
plot_heatmap(estimated_mat, legend="observed", has.legend=F,
             xnames = "Gene", ynames = "Sample",
             min=0, max=1)
par(mar = rep(2, 4))
plot(output_5_10$cv_loss, xlab="", ylab="")
plot(output_5_10$AICs, xlab="", ylab="")
