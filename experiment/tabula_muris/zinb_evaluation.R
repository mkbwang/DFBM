rm(list=ls())
library(ggplot2)
library(reshape2)

real_zeroprob <- read.table("data/tabula_muris/template/pi_hat.tsv",
                            sep='\t')

real_zeroprob_vec <- as.vector(as.matrix(real_zeroprob))

estimated_zeroprob <- read.table("data/tabula_muris/template/denoised_results/zero_inflation_probs_sign_series.tsv",
                                 sep='\t', header=T)
rownames(estimated_zeroprob) <- estimated_zeroprob$X
estimated_zeroprob$X <- NULL
estimated_zeroprob_vec <- as.vector(as.matrix(estimated_zeroprob))

zeroprob_df <- data.frame(real = real_zeroprob_vec,
                          estimate = estimated_zeroprob_vec)
zeroprob_df$difference <- zeroprob_df$estimate - zeroprob_df$real


subset_zeroprob_df <- zeroprob_df[sample(nrow(zeroprob_df), 50000), ]

library(ggplot2)
ggplot(subset_zeroprob_df, aes(x=real, y=estimate)) + geom_point(size=0.1) + 
  xlim(0, 1) + ylim(0, 1) + geom_abline(color="blue")+ xlab("Real Zero Probability") + ylab("Estimated Zero Probability")


ggplot(zeroprob_df, aes(x=difference)) + 
  geom_histogram(bins = 40, fill = "gray", color = "black", alpha = 0.7) +  # You can adjust the number of bins here
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +  # Add a vertical line at x = 0
  labs(x = "Zero Probability (Estimate - Real)",
       y = "Frequency")

