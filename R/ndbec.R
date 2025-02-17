


ndbec <- function(count_mat, quantiles=seq(0.1, 0.9, 0.1),
                  increment=0.9, max_K=10, lambdas=c(0.01, 0.1, 1)){

  prevalences <- colMeans(count_mat > 0)

  # select quantiles that are candidates for the thresholds
  unique_quantiles <- apply(count_mat, 2, quantile, probs=quantiles) |> as.vector() |>
    unique()
  unique_quantiles <- round(unique_quantiles) |> unique() |> sort()
  # browser()
  # count the number of entries larger than all the candidate thresholds
  sum_series <- rep(0, length(unique_quantiles))
  for (i in 1:length(unique_quantiles)){
    sum_series[i] <- sum(count_mat > unique_quantiles[i])
  }

  # select thresholds so that between adjacent thresholds, at least certain proportion of binary entries changed
  selected_order <- c(1)
  j <- 1
  while (j < length(unique_quantiles)){
    available_indices <- which(sum_series < sum_series[j] * increment)
    if (length(available_indices) == 0){
      break
    }
    next_j <- min(available_indices)
    selected_order <- c(selected_order, next_j)
    j <- next_j
  }
  selected_thresholds <- unique_quantiles[selected_order]
  print(sprintf("%d thresholds selected", length(selected_thresholds)))
  # browser()
  # fit logisticCF to each binary slices
  prob_mats <- list()
  ranks <- rep(0, length(selected_thresholds))

  for (j in 0:(length(selected_thresholds)-1)){

    print(sprintf("Fitting binary matrix between cutoff %d and %d", j, j+1))
    if (j == 0){
      mask <- matrix(1, nrow=nrow(count_mat), ncol=ncol(count_mat))
    } else{
      mask <- 1*(count_mat > selected_thresholds[j])
    }
    observation <- 1*(count_mat > selected_thresholds[j+1])

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
                              max_K=max_K, lambdas=lambdas)


    ranks[j+1] <- lcf_fit$selected_K
    estimated_pi <- lcf_fit$pi

    if (any(row_sums == 0) | any(col_sums == 0)){
      pi_mat[rows_retain, cols_retain] <- estimated_pi
    } else{
      pi_mat <- estimated_pi
    }

    prob_mats[[j+1]] <- pi_mat
  }

  interval_vals <- rep(0, length(selected_thresholds))
  for (j in 1:(length(selected_thresholds)-1)){
    lower_bound <- selected_thresholds[j]
    upper_bound <- selected_thresholds[j+1]

    interval_vals[j] <- mean(count_mat[count_mat > lower_bound &
                                      count_mat <= upper_bound])
  }

  # for the values larger than the largest threshold, I take median
  upper_bound <- selected_thresholds[length(selected_thresholds)]
  interval_vals[length(selected_thresholds)] <- median(count_mat[count_mat > upper_bound])

  # calculate expected counts
  expected_counts <- matrix(0, nrow=nrow(count_mat), ncol=ncol(count_mat))
  for (j in 1:length(selected_thresholds)){

    if (j==1){
      sprob1 <- prob_mats[[1]]
    } else{
      sprob1 <- sprob1*prob_mats[[j]]
    }

    if (j==length(selected_thresholds)){
      sprob2 <- 0
    } else{
      sprob2 <- sprob1*prob_mats[[j+1]]
    }
    expected_counts <- expected_counts + interval_vals[j] * (sprob1 - sprob2)

  }

  return(list(thresholds=selected_thresholds,
              optimal_ranks = ranks,
              prob_mats = prob_mats,
              denoised_counts=expected_counts))


}

