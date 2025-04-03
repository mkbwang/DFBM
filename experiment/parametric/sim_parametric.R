
folder <- "experiment/parametric/"

counts_mat <- read.csv(file.path(folder,"nb_counts.csv"),
                       row.names=1) |> as.matrix()



dfbm_denoise <- dfbm(counts_mat, increment=0.7,
                     max_K=8, ncores=4)

dfbm_denoise_count_df <- data.frame(dfbm_denoise$denoised_counts)

write.csv(dfbm_denoise_count_df, "experiment/parametric/nb_denoised_dfbm.csv",
          row.names=F)

