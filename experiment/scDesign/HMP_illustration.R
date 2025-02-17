

rm(list=ls())
library(reshape2)

# load the true distribution parameters
parameters <- readRDS("experiment/scDesign/HMP/fitted_parameters.rds")
zinb_parameters <- parameters$ZINB

plaque_mean <- rep(0, 500)
plaque_mean[setdiff(1:500, zinb_parameters$Plaque$gene_sel3)] <-
  zinb_parameters$Plaque$marginal_param1[, 3]

tongue_mean <- rep(0, 500)
tongue_mean[setdiff(1:500, zinb_parameters$Tongue$gene_sel3)] <-
  zinb_parameters$Tongue$marginal_param1[, 3]

throat_mean <- rep(0, 500)
throat_mean[setdiff(1:500, zinb_parameters$Throat$gene_sel3)] <-
  zinb_parameters$Throat$marginal_param1[, 3]

# find the top 100 genes with the largest variance
zinb_means <- cbind(plaque_mean, tongue_mean, throat_mean)
meanval_var <-apply(zinb_means,1, var)
var_ranking <- sort(meanval_var, index.return=TRUE, decreasing=TRUE)

selected_columns <- var_ranking$ix[1:100]



zinb_counts <- readRDS("experiment/scDesign/HMP/ZINB_sim.rds") |> t()

metadata_zinb <- read.csv("experiment/scDesign/HMP/metadata_zinb.csv",
                            row.names=1)

plaque_samples <- sample(which(metadata_zinb$Type == "Plaque"), 30)
tongue_samples <- sample(which(metadata_zinb$Type == "Tongue"), 30)
throat_samples <- sample(which(metadata_zinb$Type == "Throat"), 30)
selected_samples <- c(plaque_samples, tongue_samples, throat_samples)


zinb_count_subset <- zinb_counts[selected_samples, selected_columns]

plot_heatmap <- function(input_matrix, mask = NULL, color.low="black", color.high="white",
                         color.na="#166919", min=0, max=1){
  rownames(input_matrix) <- NULL
  colnames(input_matrix) <- NULL

  if(!is.null(mask)){
    rownames(mask) <- NULL
    colnames(mask) <- NULL
    input_matrix[mask == 0] <- NA
  }

  long_matrix <- melt(input_matrix)

  plot_object <- ggplot(data = long_matrix, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() + labs(x = NULL, y=NULL, fill=NULL) +
    scale_fill_gradient2(low=color.low, high=color.high, na.value=color.na,
                         limits=c(min, max), oob=scales::squish)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "null"),
          panel.spacing = unit(c(0, 0, 0, 0), "null"))

  return(plot_object)

}

original_plot <- plot_heatmap(zinb_count_subset,color.low="white", color.high="#db1456",
                              min=0, max=100)


zinb_denoise_4 <- ndbec(count_mat=zinb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.4,
                        max_K=10,
                        lambdas=0.1)

denoise_mat <- zinb_denoise_4$denoised_counts
denoise_count_subset <- denoise_mat[selected_samples, selected_columns]

denoised_plot <- plot_heatmap(denoise_count_subset,color.low="white", color.high="#db1456",
                              min=0, max=100)

threshold_values <- zinb_denoise_4$thresholds
binary_slices <- list()
for (j in 1:length(threshold_values)){
  binary_slices[[j]] <-1*( zinb_count_subset > threshold_values[j])
}
binary_slice_plot <- list()
for (j in 1:length(threshold_values)){
  binary_slice_plot[[j]] <- plot_heatmap(binary_slices[[j]],color.low="white", color.high="#5c5e5d",
                                         min=0, max=1)
}

binary_slice_withNA_plot <- list()
for (j in 1:(length(threshold_values)-1)){
  binary_slice_withNA_plot[[j]] <- plot_heatmap(binary_slices[[j+1]], mask = binary_slices[[j]],
                                           color.low="white", color.high="#5c5e5d", color.na="#4648c2",
                                           min=0, max=1)
}


estim_probs <- list()
for (j in 1:length(threshold_values)){
  if(j == 1){
    estim_probs[[j]] <- zinb_denoise_4$prob_mats[[1]]
  } else{
    estim_probs[[j]] <- estim_probs[[j-1]] * zinb_denoise_4$prob_mats[[j]]
  }
}
for (j in 1:length(threshold_values)){
  estim_probs[[j]] <- estim_probs[[j]][selected_samples, selected_columns]
}

estim_probs_plot <- list()
for (j in 1:length(threshold_values)){
  estim_probs_plot[[j]] <- plot_heatmap(estim_probs[[j]],color.low="white", color.high="#1a24d9",
                                         min=0, max=1)
}


thresholds <- zinb_denoise_4$thresholds


prob_example_1 <- c(estim_probs[[1]][1,3], estim_probs[[2]][1,3],
                    estim_probs[[3]][1,3], estim_probs[[4]][1,3],
                    estim_probs[[5]][1,3], estim_probs[[6]][1,3])

prob_example_2 <- c(estim_probs[[1]][6,1], estim_probs[[2]][6,1],
                    estim_probs[[3]][6,1], estim_probs[[4]][6,1],
                    estim_probs[[5]][6,1], estim_probs[[6]][6,1])


prob_example_3 <- c(estim_probs[[1]][90,1], estim_probs[[2]][90,1],
                    estim_probs[[3]][90,1], estim_probs[[4]][90,1],
                    estim_probs[[5]][90,1], estim_probs[[6]][90,1])


plot_dist <- function(thresholds, probs){
  df <- data.frame(Thresholds = c(thresholds[1], rep(thresholds[2:6], each=2), thresholds[6] + 10),
                             Probs = rep(probs, each=2))

  generated_plot <- ggplot(df, aes(x=Thresholds, y=Probs)) +
    geom_line() +
    labs(x = NULL, y=NULL, fill=NULL) +
    xlim(0, thresholds[6] + 10)+ ylim(0, 1)+
    theme_classic()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "null"),
          panel.spacing = unit(c(0, 0, 0, 0), "null"))

  return(generated_plot)
}


prob_plot_1 <- plot_dist(thresholds, prob_example_1)

prob_plot_2 <- plot_dist(thresholds, prob_example_2)

prob_plot_3 <- plot_dist(thresholds, prob_example_3)


