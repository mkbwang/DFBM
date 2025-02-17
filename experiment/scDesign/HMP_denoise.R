
rm(list=ls())

# Load the data

nb_counts <- readRDS("experiment/scDesign/HMP/NB_sim.rds") |> t()
metadata_nb <- data.frame(Sample=sprintf("Sample%d", seq(1, nrow(nb_counts))),
                          Type=rownames(nb_counts))
rownames(nb_counts) <- metadata_nb$Sample
write.table(nb_counts, "experiment/scDesign/HMP/NB_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)
write.csv(metadata_nb, "experiment/scDesign/HMP/metadata_nb.csv",
          row.names=FALSE, quote=FALSE)

poisson_counts <- readRDS("experiment/scDesign/HMP/poisson_sim.rds") |> t()
metadata_poisson <- data.frame(Sample=sprintf("Sample%d", seq(1, nrow(poisson_counts))),
                               Type=rownames(poisson_counts))
rownames(poisson_counts) <- metadata_poisson$Sample
write.table(poisson_counts, "experiment/scDesign/HMP/poisson_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)
write.csv(metadata_poisson, "experiment/scDesign/HMP/metadata_poisson.csv",
          row.names=FALSE, quote=FALSE)


zinb_counts <- readRDS("experiment/scDesign/HMP/ZINB_sim.rds") |> t()
metadata_zinb <- data.frame(Sample=sprintf("Sample%d", seq(1, nrow(zinb_counts))),
                            Type=rownames(zinb_counts))
rownames(zinb_counts) <- metadata_zinb$Sample
write.table(zinb_counts, "experiment/scDesign/HMP/zinb_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)
write.csv(metadata_zinb, "experiment/scDesign/HMP/metadata_zinb.csv",
          row.names=FALSE, quote=FALSE)


full_counts <- readRDS("experiment/scDesign/HMP/ZINB_sim.rds") |> t()
write.table(full_counts, "experiment/scDesign/HMP/full_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)

nb_denoise_1 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.1,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_1$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_1.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_2 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.2,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_2$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_2.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_3 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.3,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_3$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_3.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_4 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.4,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_4$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_4.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_5 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.5,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_5$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_5.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_6 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.6,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_6$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_6.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_7 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.7,
                      max_K=10,
                      lambdas=0.1)

write.table(nb_denoise_7$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_7.txt",
            sep='\t', quote=FALSE, row.names=TRUE)

nb_denoise_8 <- ndbec(count_mat=nb_counts,
                    quantiles=seq(0.1, 0.9, 0.1),
                    increment=0.8,
                    max_K=10,
                    lambdas=0.1)

write.table(nb_denoise_8$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_8.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


nb_denoise_9 <- ndbec(count_mat=nb_counts,
                    quantiles=seq(0.1, 0.9, 0.1),
                    increment=0.9,
                    max_K=10,
                    lambdas=0.1)

write.table(nb_denoise_9$denoised_counts, file="experiment/scDesign/HMP/NDBEC/nb_denoise_9.txt",
            sep='\t', quote=FALSE, row.names=TRUE)





zinb_denoise_1 <- ndbec(count_mat=nb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.1,
                        max_K=10,
                        lambdas=0.1)

write.table(zinb_denoise_1$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_1.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


zinb_denoise_2 <- ndbec(count_mat=nb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.2,
                      max_K=10,
                      lambdas=0.1)

write.table(zinb_denoise_2$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_2.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


zinb_denoise_3 <- ndbec(count_mat=zinb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.3,
                        max_K=10,
                        lambdas=0.1)

write.table(zinb_denoise_3$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_3.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


zinb_denoise_4 <- ndbec(count_mat=zinb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.4,
                        max_K=10,
                        lambdas=0.1)


write.table(zinb_denoise_4$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_4.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


zinb_denoise_5 <- ndbec(count_mat=zinb_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.5,
                      max_K=10,
                      lambdas=0.1)

write.table(zinb_denoise_5$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_5.txt",
            sep='\t', quote=FALSE, row.names=TRUE)



zinb_denoise_6 <- ndbec(count_mat=zinb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.6,
                        max_K=10,
                        lambdas=0.1)

write.table(zinb_denoise_6$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_6.txt",
            sep='\t', quote=FALSE, row.names=TRUE)



zinb_denoise_7 <- ndbec(count_mat=zinb_counts,
                         quantiles=seq(0.1, 0.9, 0.1),
                         increment=0.7,
                         max_K=10,
                         lambdas=0.1)

write.table(zinb_denoise_7$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_7.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


zinb_denoise_8 <- ndbec(count_mat=zinb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.8,
                        max_K=10,
                        lambdas=0.1)

write.table(zinb_denoise_8$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_8.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


zinb_denoise_9 <- ndbec(count_mat=zinb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.9,
                        max_K=10,
                        lambdas=0.1)

write.table(zinb_denoise_9$denoised_counts, file="experiment/scDesign/HMP/NDBEC/zinb_denoise_9.txt",
            sep='\t', quote=FALSE, row.names=TRUE)



# start <- proc.time()
# full_denoise <- ndbec(count_mat=full_counts,
#                       quantiles=seq(0.1, 0.9, 0.1),
#                       increment=0.8,
#                       max_K=10,
#                       lambdas=0.1)
# end <- proc.time()
# write.table(full_denoise$denoised_counts, file="experiment/scDesign/HMP/full_denoise.txt",
#             sep='\t', quote=FALSE, row.names=TRUE)




