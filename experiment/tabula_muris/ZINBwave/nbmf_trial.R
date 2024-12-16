
rm(list=ls())
count_mat <- read.table("experiment/tabula_muris/ZINBwave/simulated_counts.tsv",
                        header=TRUE, row.names=1, sep="\t")

library(Matrix)

example_slice <- 1*(count_mat > 0)
example_mask <- matrix(1, nrow=nrow(example_slice), ncol=ncol(example_slice))

## TODO: visualize the entries at risk and the binary result

cv_result <- cv.nbmf(Y = example_slice, Z = example_mask, max_k=10)
AICs <- rep(0, 10)
for (j in 1:10){
  print(j)
  result <- nbmf(Y = example_slice, Z = example_mask, k=j)
  end <- proc.time()
  AICs[j] <- result$AIC
}

# TODO: plot AICs and cv loss

selected_result <- nbmf(Y = example_slice, Z = example_mask,
                        k=cv_result$best_rank)


# TODO: visualize the estimated probability

output <- list(cv_loss = cv_result$loss, AICs=AICs,
               fit = selected_result)

saveRDS(output, "experiment/tabula_muris/ZINBwave/nbmf_zeroprob_output.rds")




