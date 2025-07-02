

#' Denoise count matrices by factorization of binary masks
#'
#' @param count_mat count matrix, nsamples * nfeatures
#' @param quantiles the candidate thresholds are selected from the quantiles of each feature values
#' @param increment a value between 0 and 1, at least certain percentage of entries larger than one threshold is smaller than the next threshold
#' @param cutoffs if a vector of cutoffs are provided, then no need to derive thresholds
#' @param fix_Ks the user can fix the rank number for all the binary masks
#' @param max_K maximum number of ranks for binary matrix factorization
#' @param lambdas ridge penalty parameters to try for the entries in the factorized matrices
#' @param ignore if less than a certain proportion of entries are "observed" for a column, that column is ignored for matrix factorization
#' @param cap only denoise entries whose values are smaller than this value, by default 99% quantile of all entries if NA provided
#' @param ncores number of cores for parallel computing, default 1
#'
#' @importFrom stats quantile median
#' @returns a list object that includes the denoised counts, denoised probability matrices,thresholds and optimal ranks
#' @export
dfbm <- function(count_mat, quantiles=seq(0.1, 0.9, 0.1),
                  increment=0.9, cutoffs=NULL, fix_Ks=NULL, max_K=10, lambdas=c(0.01, 0.1, 1),
                 ignore=0, cap=NA,
                 ncores=1){

  nsample <- nrow(count_mat)
  nfeature <- ncol(count_mat)
  maximum_threshold <- quantile(as.vector(count_mat), 0.99)
  if (!is.null(cap)){
    maximum_threshold <- cap
  }

  prevalences <- colMeans(count_mat > 0)

  if (is.null(cutoffs)){
    # select quantiles that are candidates for the thresholds
    unique_quantiles <- apply(count_mat, 2, quantile, probs=quantiles) |> as.vector() |>
      unique() |> sort()
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
  } else{
    selected_thresholds <- cutoffs
  }
  # remove automatically decided thresholds that were larger than the maximum value
  selected_thresholds <- selected_thresholds[selected_thresholds <= maximum_threshold]

  if (!is.null(fix_Ks)){
    fix_Ks <- fix_Ks[1:length(selected_thresholds)]
  }

  print(sprintf("%d thresholds selected", length(selected_thresholds)))

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
    cols_retain <- seq(1, ncol(count_mat))

    # Do I need to subset rows and columns
    col_means <- colMeans(mask)
    if (any(col_means <= ignore)){ # remove genes which are already very censored
      cols_retain <- which(col_means > ignore)
    }
    mask <- mask[, cols_retain, drop=FALSE]
    row_means <- rowMeans(mask)
    if (any(row_means == 0)){ # remove rows are already all censored
      rows_retain <- which(row_means > 0)
    }
    mask <- mask[rows_retain, , drop=FALSE]
    observation <- observation[rows_retain, cols_retain, drop=FALSE]

    if (nrow(observation) == 1 | ncol(observation) == 1){ # when the dimension of matrix reduce to 1, we can only take the mean
      pi_scalar <- sum(observation) / sum(mask)
      estimated_pi <- matrix(pi_scalar, nrow=nrow(observation), ncol=ncol(observation))
    } else{
      # logistic collaborative filtering for denoising binary matrices
      if (is.null(fix_Ks)){ # need to select rank
        lcf_fit <- cv.logisticcfR(X=observation, Z=mask,
                                  max_K=max_K, lambdas=lambdas, ncores=ncores)
        ranks[j+1] <- lcf_fit$selected_K
        estimated_pi <- lcf_fit$pi
      } else{
        lcf_fit <- logisticcfR(X=observation, Z=mask, K=fix_Ks[j+1], lambda=sample(lambdas, 1))
        ranks[j+1] <- fix_Ks[j+1]
        estimated_pi <- lcf_fit$pi
      }
    }

    if (any(row_means <= ignore) | any(col_means <= ignore)){
      pi_mat[rows_retain, cols_retain] <- estimated_pi
    } else{
      pi_mat <- estimated_pi
    }
    prob_mats[[j+1]] <- pi_mat
  }

  cap_vals <- apply(count_mat, 2, quantile, probs=1-ignore/2)
  cap_vals_mat <- matrix(rep(cap_vals, nsample), ncol=nsample) |> t()

  interval_vals <- rep(0, length(selected_thresholds))
  for (j in 1:(length(selected_thresholds)-1)){
    lower_bound <- selected_thresholds[j]
    upper_bound <- selected_thresholds[j+1]

    interval_vals[j] <- mean(count_mat[count_mat > lower_bound &
                                      count_mat <= upper_bound])
  }


  upper_bound <- selected_thresholds[length(selected_thresholds)]
  interval_vals[length(selected_thresholds)] <- mean(count_mat[count_mat > upper_bound &
                                                                 count_mat <= maximum_threshold])

  # calculate expected counts
  expected_counts <- matrix(0, nrow=nrow(count_mat), ncol=ncol(count_mat))
  for (j in 1:length(selected_thresholds)){

    prev_threshold <- selected_thresholds[j]
    if (j==1){
      sprob1 <- prob_mats[[1]]
    } else{
      sprob1 <- sprob1*prob_mats[[j]]
    }

    if (j==length(selected_thresholds)){
      sprob2 <- matrix(0, nrow=nsample, ncol=nfeature)
    } else{
      sprob2 <- sprob1*prob_mats[[j+1]]
    }

    # if (interpolate){
    #   tail_values_mat <- (cap_vals_mat + prev_threshold) / 2
    #   estim_tail_mask <- sprob2 == 0 # whether probability of larger than the next threshold is zero
    #   if(any(estim_tail_mask)){ # some entries' probability drop to zero
    #     expected_counts[estim_tail_mask] <- expected_counts[estim_tail_mask] +
    #       tail_values_mat[estim_tail_mask] * (sprob1[estim_tail_mask] - sprob2[estim_tail_mask])
    #   }
    #   if(!all(estim_tail_mask)){
    #     expected_counts[!estim_tail_mask] <- expected_counts[!estim_tail_mask] +
    #       interval_vals[j] * (sprob1[!estim_tail_mask] - sprob2[!estim_tail_mask])
    #   }
    # } else{
    expected_counts <- expected_counts + interval_vals[j] * (sprob1 - sprob2)
    # }
  }

  # do not change entries with large values
  expected_counts[count_mat > maximum_threshold] <- count_mat[count_mat > maximum_threshold]

  return(list(thresholds=selected_thresholds,
              optimal_ranks = ranks,
              prob_mats = prob_mats,
              denoised_counts=expected_counts))


}

