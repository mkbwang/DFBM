
rm(list=ls())
library(phyloseq)
data <- readRDS("experiment/COHRA2/phyasv_visit24.rds")

metadata <- sample_data(data) |> as.data.frame()
sample_filter <- grepl("pre-incident", metadata$CaseStatus)
metadata_subset <- metadata[sample_filter, ]
counts <- otu_table(data)@.Data
counts_subset <- counts[sample_filter, ]

# metadata_yr1 <- metadata[metadata$Visit == 5, ]
# counts <- counts[metadata$Visit == 5, ]
prevalences <- colMeans(counts_subset > 0)
counts_subset <- counts_subset[, prevalences > 0.05]

# denoised_output <- ndbec(count_mat=counts_subset,
#                          quantiles=seq(0.1, 0.9, 0.1),
#                          increment=0.6,
#                          max_K=10,
#                          lambdas=0.1)

thresholds <- rep(0, 9)
for (j in 1:9){
  print(j)
  denoised_output <- ndbec(count_mat=counts_subset,
                           quantiles=seq(0.1, 0.9, 0.1),
                           increment=j*0.1,
                           max_K=10,
                           lambdas=0.1)
  thresholds[j] <- length(denoised_output$thresholds)
  write.csv(as.data.frame(denoised_output$denoised_counts),
            sprintf("experiment/COHRA2/denoise/16S_yr2_%d.csv", j),
            quote=F, row.names = F)

}




# export metadata and expected counts

write.csv(data.frame(metadata_subset),
          "experiment/COHRA2/metadata_yr2.csv",
          quote=F)


write.csv(as.data.frame(counts_subset),
          "experiment/COHRA2/16S_counts_yr2.csv",
          quote=F)


