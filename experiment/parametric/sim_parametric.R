
rm(list=ls())

folder <- "experiment/parametric/"

counts_mat <- read.csv(file.path(folder,"nb_count_10.csv"),
                       row.names=1) |> as.matrix()

thresholds <- c(0,4,8,12,16,20,25,31,38,45,52,
                58,64,70,76,83,89,96,103,109,117,129,148)


dfbm_denoise <- dfbm(counts_mat, cutoffs=thresholds, ignore=0,
                     fix_Ks=rep(3, length(thresholds)), ncores=4)

denoised_counts <- dfbm_denoise$denoised_counts


dfbm_denoise_1 <- dfbm(counts_mat, cutoffs=thresholds, ignore=0.1,
                       fix_Ks=rep(3, length(thresholds)), max_K=8, ncores=4)

denoised_counts_1 <- dfbm_denoise_1$denoised_counts



