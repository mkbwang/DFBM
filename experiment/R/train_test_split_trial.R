
RhpcBLASctl::blas_set_num_threads(1)
library(ggplot2)
# load all the compositions
# flocyt_composition <- read.csv("inst/extdata/flocyt_proj_composition_leaf.csv")
flocyt_composition <- read.csv("inst/extdata/flocyt_proj_composition_leaf_merged.csv")


begin <- proc.time()
cv_dfbm_result <- dfbm(A=flocyt_composition,
                              alpha_grid=seq(1, 0.1, -0.1), ncores=8, n_splits=5,
                              tol=2e-5)
cv_dfbm_denoised <- cv_dfbm_result$denoised |> as.data.frame()
cv_time <- proc.time() - begin


# load CV diagnostics
cv_diagnostics <- cv_dfbm_result$cv
RPS_plot <- ggplot(cv_diagnostics, aes(x = alpha, y = rps)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = rps - rps_se, ymax = rps + rps_se),
                width = 0.03, colour = "black") +
  geom_point(size = 1.5, colour = "black") +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  labs(x = "Alpha", y = "RPS", title = "RPS with standard error") +
  theme_bw()

CRPS_plot <- ggplot(cv_diagnostics, aes(x = alpha, y = crps)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = crps - crps_se, ymax = crps + crps_se),
                width = 0.03, colour = "black") +
  geom_point(size = 1.5, colour = "black") +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  labs(x = "Alpha", y = "CRPS", title = "CRPS with standard error") +
  theme_bw()

MSE_plot <- ggplot(cv_diagnostics, aes(x = alpha, y = mse)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = mse - mse_se, ymax = mse + mse_se),
                width = 0.03, colour = "black") +
  geom_point(size = 1.5, colour = "black") +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  labs(x = "Alpha", y = "MSE", title = "MSE with standard error") +
  theme_bw()


# pick alpha=0.5 and give a final run
dfbm_final <- dfbm(A=flocyt_composition, alpha=0.5)

dfbm_final_denoised <- dfbm_final$denoised |> as.data.frame()

write.csv(dfbm_final_denoised, "experiment/data/denoised_composition_leaf_merged.csv", row.names=FALSE)

cd4naive <- cbind(dfbm_final_denoised$pcd4n_adj_pct, flocyt_composition$pcd4n_adj_pct) |> as.data.frame()
cd4temra <- cbind(dfbm_final_denoised$pcd4temra_adj_pct, flocyt_composition$pcd4temra_adj_pct) |> as.data.frame()
cd8naive <- cbind(dfbm_final_denoised$pcd8n_adj_pct, flocyt_composition$pcd8n_adj_pct) |> as.data.frame()
cd8temra <- cbind(dfbm_final_denoised$pcd8temra_adj_pct, flocyt_composition$pcd8temra_adj_pct) |> as.data.frame()
naiveb <- cbind(dfbm_final_denoised$pnaiveb_adj_pct, flocyt_composition$pnaiveb_adj_pct) |> as.data.frame()
igdminus <- cbind(dfbm_final_denoised$pigd_minus_memb_adj_pct, flocyt_composition$pigd_minus_memb_adj_pct) |> as.data.frame()
igdplus <- cbind(dfbm_final_denoised$pigd_plus_memb_adj_pct, flocyt_composition$pigd_plus_memb_adj_pct) |> as.data.frame()

colnames(cd4naive) <- colnames(cd4temra) <- colnames(cd8naive) <-
  colnames(cd8temra) <- colnames(naiveb) <- colnames(igdminus) <-
  colnames(igdplus) <- c("Denoised", "Original")

comparison_single_plot <- function(df, title){
  max_value <- max(df$Original)
  pl <- ggplot(df, aes(x = Original, y = Denoised)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "red") +
    labs(x = "Original", y = "Denoised", title=title) +
    scale_x_continuous(limits = c(0, max_value)) +
    scale_y_continuous(limits = c(0, max_value)) +
    theme_bw()

  pl
}
cd4naive_plot <- comparison_single_plot(df=cd4naive, title="CD4 Naive T Cell")
cd4temra_plot <- comparison_single_plot(df=cd4temra, title="CD4 TemRA Cell")
cd8naive_plot <- comparison_single_plot(df=cd8naive, title="CD8 Naive T Cell")
cd8temra_plot <- comparison_single_plot(df=cd8naive, title="CD8 TemRA Cell")
naiveb_plot <- comparison_single_plot(df=naiveb, title="Naive B Cell")
igdminus_plot <- comparison_single_plot(df=igdminus, title="IgD- B Cell")
igdplus_plot <- comparison_single_plot(df=igdminus, title="IgD+ B Cell")

max_value <- max(cd4naive$Original)



# begin <- proc.time()
# dfbm_single_result_temp <- dfbm(A=flocyt_composition,alpha=0.8)
# time5 <- proc.time() - begin

#
# write.csv(dfbm_denoised, "experiment/data/denoised_composition_leaf_merged.csv", row.names=FALSE)


#TODO: the rest are obsolete for now
# split into train and test
# for (j in 1:20){
#   print(j)
#   set.seed(j)
#   train_ids <- sample(seq(1:nrow(flocyt_composition)),
#                       0.7*nrow(flocyt_composition))
#   test_ids <- setdiff(seq(1:nrow(flocyt_composition)), train_ids)
#
#   train_compositions <- flocyt_composition[train_ids, ]
#   test_compositions <- flocyt_composition[test_ids, ]
#
#
#   begin <- proc.time()
#   train_dfbm_result <- dfbm(A=train_compositions, alpha=0.55)
#   end <- proc.time()
#   train_dfbm_denoised <- train_dfbm_result$denoised |> as.data.frame()
#
#
#   test_dfbm_result <- predict.dfbm(train_dfbm_result, newdata=test_compositions)
#   test_dfbm_denoised <- test_dfbm_result$denoised |> as.data.frame()
#
#
#   denoise_result <- list(train_ids=train_ids, test_ids=test_ids,
#                          train_original=train_compositions,
#                          test_original=test_compositions,
#                          train_denoised=train_dfbm_denoised,
#                          test_denoised=test_dfbm_denoised)
#   saveRDS(denoise_result, sprintf("experiment/data/denoise_flocyt_merged_%d.rds", j))
#
# }
#
#
#
#
# flocyt_composition <- read.csv("inst/extdata/flocyt_proj_composition_leaf.csv")
# for (j in 1:20){
#   print(j)
#   set.seed(j)
#   train_ids <- sample(seq(1:nrow(flocyt_composition)),
#                       0.7*nrow(flocyt_composition))
#   test_ids <- setdiff(seq(1:nrow(flocyt_composition)), train_ids)
#
#   train_compositions <- flocyt_composition[train_ids, ]
#   test_compositions <- flocyt_composition[test_ids, ]
#
#
#   begin <- proc.time()
#   train_dfbm_result <- dfbm(A=train_compositions, alpha=0.55)
#   end <- proc.time()
#   train_dfbm_denoised <- train_dfbm_result$denoised |> as.data.frame()
#
#
#   test_dfbm_result <- predict.dfbm(train_dfbm_result, newdata=test_compositions)
#   test_dfbm_denoised <- test_dfbm_result$denoised |> as.data.frame()
#
#
#   denoise_result <- list(train_ids=train_ids, test_ids=test_ids,
#                          train_original=train_compositions,
#                          test_original=test_compositions,
#                          train_denoised=train_dfbm_denoised,
#                          test_denoised=test_dfbm_denoised)
#   saveRDS(denoise_result, sprintf("experiment/data/denoise_flocyt_unmerged_%d.rds", j))
#
# }
#
#
#
