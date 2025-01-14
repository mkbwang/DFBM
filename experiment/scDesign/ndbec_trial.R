
# load the count matrix
count_mat <- read.table("experiment/scDesign/10xgenomics_simulated_count.tsv",
                        header=TRUE, row.names=1, sep="\t") |> as.matrix()

# calculate quantiles for each column
quantiles <- apply(count_mat, 2, quantile, probs=seq(0.1, 0.9, 0.1)) |>
  as.vector() |> unique()
quantiles <- round(quantiles) |> unique() |> sort()

sum_series <- rep(0, length(quantiles))
for (i in 1:length(quantiles)){
  sum_series[i] <- sum(count_mat > quantiles[i])
}

selected_order <- c(1)
j <- 1
while (j < length(quantiles)){
  next_j <- min(which(sum_series < sum_series[j] * 0.8))
  selected_order <- c(selected_order, next_j)
  j <- next_j
}
selected_quantiles <- quantiles[selected_order]


# begin fitting logisticcf to each slice

prob_mats <- list()
ranks <- rep(0, length(selected_quantiles))

begin <- proc.time()
for (j in 0:(length(selected_quantiles)-1)){

  print(sprintf("Fitting binary matrix between cutoff %d and %d", j, j+1))
  if (j == 0){
    mask <- matrix(1, nrow=nrow(count_mat), ncol=ncol(count_mat))
  } else{
    mask <- 1*(count_mat > selected_quantiles[j])
  }
  observation <- 1*(count_mat > selected_quantiles[j+1])

  pi_mat <- matrix(0, nrow=nrow(count_mat), ncol=ncol(count_mat))
  rows_retain <- seq(1, nrow(count_mat))
  ncols_retain <- seq(1, ncol(count_mat))

  # Do I need to subset rows and columns
  row_sums <- rowSums(mask)
  if (any(row_sums == 0)){
    rows_retain <- which(row_sums > 0)
    mask <- mask[rows_retain, ]
    observation <- observation[rows_retain, ]
  }

  col_sums <- colSums(mask)
  if (any(col_sums == 0)){
    cols_retain <- which(col_sums > 0)
    mask <- mask[, cols_retain]
    observation <- observation[, cols_retain]
  }

  lcf_fit <- cv.logisticcfR(X=observation, Z=mask,
                            max_K=10, lambdas=0.1)

  ranks[j+1] <- lcf_fit$selected_K
  estimated_pi <- lcf_fit$pi

  if (any(row_sums == 0) | any(col_sums == 0)){
    pi_mat[rows_retain, cols_retain] <- estimated_pi
  } else{
    pi_mat <- estimated_pi
  }

  prob_mats[[j+1]] <- pi_mat

}
end <- proc.time() - begin

# save the results
saveRDS(prob_mats, "experiment/scDesign/ndbec_trial_prob_mats.rds")
interval_vals <- rep(0, length(selected_quantiles))


for (j in 1:(length(selected_quantiles)-1)){
  lower_bound <- selected_quantiles[j]
  upper_bound <- selected_quantiles[j+1]

  interval_vals[j] <- mean(count_mat[count_mat > lower_bound &
                                       count_mat <= upper_bound])
}
# for the values larger than the largest threshold, I take median
upper_bound <- selected_quantiles[length(selected_quantiles)]
interval_vals[length(selected_quantiles)] <- median(count_mat[count_mat > upper_bound])

# calculate expected counts
expected_counts <- matrix(0, nrow=nrow(count_mat), ncol=ncol(count_mat))
for (j in 1:length(selected_quantiles)){

  if (j==1){
    sprob1 <- prob_mats[[1]]
  } else{
    sprob1 <- sprob1*prob_mats[[j]]
  }

  if (j==length(selected_quantiles)){
    sprob2 <- 0
  } else{
    sprob2 <- sprob1*prob_mats[[j+1]]
  }
  expected_counts <- expected_counts + interval_vals[j] * (sprob1 - sprob2)

}


prob_gt0 <- prob_mats[[1]]

write.table(data.frame(prob_gt0),
            "experiment/scDesign/ndbec_trial_prob_gt0.tsv",
            sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

write.table(data.frame(expected_counts),
            "experiment/scDesign/ndbec_trial_expected_counts.tsv",
            sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)


