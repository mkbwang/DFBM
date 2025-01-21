
rm(list=ls())

# Load the data

nb_counts <- readRDS("experiment/scDesign/HMP/NB_sim.rds") |> t()
write.table(nb_counts, "experiment/scDesign/HMP/NB_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)

poisson_counts <- readRDS("experiment/scDesign/HMP/poisson_sim.rds") |> t()
write.table(poisson_counts, "experiment/scDesign/HMP/poisson_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)

zinb_counts <- readRDS("experiment/scDesign/HMP/ZINB_sim.rds") |> t()
write.table(zinb_counts, "experiment/scDesign/HMP/zinb_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)

full_counts <- readRDS("experiment/scDesign/HMP/ZINB_sim.rds") |> t()
write.table(full_counts, "experiment/scDesign/HMP/full_sim.txt", sep="\t",
            quote=FALSE, row.names=TRUE)


denoise_output <- ndbec(count_mat=nb_counts,
                        quantiles=seq(0.1, 0.9, 0.1),
                        increment=0.8,
                        max_K=10,
                        lambdas=0.1)


nb_denoise <- ndbec(count_mat=nb_counts,
                    quantiles=seq(0.1, 0.9, 0.1),
                    increment=0.8,
                    max_K=10,
                    lambdas=0.1)

write.table(nb_denoise$denoised_counts, file="experiment/scDesign/HMP/nb_denoise.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


start <- proc.time()
poisson_denoise <- ndbec(count_mat=poisson_counts,
                         quantiles=seq(0.1, 0.9, 0.1),
                         increment=0.8,
                         max_K=10,
                         lambdas=0.1)
end <- proc.time()

write.table(poisson_denoise$denoised_counts, file="experiment/scDesign/HMP/poisson_denoise.txt",
            sep='\t', quote=FALSE, row.names=TRUE)


start <- proc.time()
zinb_denoise <- ndbec(count_mat=zinb_counts,
                         quantiles=seq(0.1, 0.9, 0.1),
                         increment=0.8,
                         max_K=10,
                         lambdas=0.1)
end <- proc.time()
write.table(zinb_denoise$denoised_counts, file="experiment/scDesign/HMP/zinb_denoise.txt",
            sep='\t', quote=FALSE, row.names=TRUE)



start <- proc.time()
full_denoise <- ndbec(count_mat=full_counts,
                      quantiles=seq(0.1, 0.9, 0.1),
                      increment=0.8,
                      max_K=10,
                      lambdas=0.1)
end <- proc.time()
write.table(full_denoise$denoised_counts, file="experiment/scDesign/HMP/full_denoise.txt",
            sep='\t', quote=FALSE, row.names=TRUE)




