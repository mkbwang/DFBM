rm(list=ls())

library(scDesign2)
library(HMP16SData)
library(dplyr)

metadata <- V35() %>% table_one()
metadata_oral <- metadata %>% filter(`HMP Body Site` == "Oral")

V35_oral <- V35() %>% subset(select = HMP_BODY_SUBSITE %in%
                               c("Subgingival Plaque", "Tongue Dorsum", "Throat"))

oral_metadata <- colData(V35_oral) |> data.frame()
counts_mat <- assay(V35_oral, "16SrRNA")
libsizes <- colSums(counts_mat)
sample_filter <- libsizes > 2000 & libsizes < 10000
prevalences <- rowMeans(counts_mat > 0)
taxa_filter <-  prevalences > 0.05
counts_mat <- counts_mat[taxa_filter, sample_filter]

# rank taxa by variance
count_var <- apply(counts_mat, 1, var)
variance_ranks <- rank(-count_var)
subset_otus <- which(variance_ranks <=500)

# only select the most variable OTUs
counts_mat <- counts_mat[subset_otus, ]
prevalences <- rowMeans(counts_mat > 0)
metadata <- oral_metadata[sample_filter, ]
metadata$HMP_BODY_SUBSITE[metadata$HMP_BODY_SUBSITE == "Subgingival Plaque"] <- "Plaque"
metadata$HMP_BODY_SUBSITE[metadata$HMP_BODY_SUBSITE == "Tongue Dorsum"] <- "Tongue"

# simulate from this template using scDesign
colnames(counts_mat) <- metadata$HMP_BODY_SUBSITE
cell_type_prop <- table(metadata$HMP_BODY_SUBSITE) / nrow(metadata)

poisson_fit <- fit_model_scDesign2(data_mat = counts_mat,
                                   cell_type_sel=c("Plaque", "Tongue", "Throat"),
                                   marginal="poisson",
                                   zp_cutoff=1,
                                   ncores=2)


poisson_sim <- simulate_count_scDesign2(model_params=poisson_fit,
                                        n_cell_new = 200,
                                        sim_method="copula",
                                        cell_type_prop = cell_type_prop)

saveRDS(poisson_sim, file="experiment/scDesign/HMP/poisson_sim.rds")


NB_fit <- fit_model_scDesign2(data_mat = counts_mat,
                              cell_type_sel=c("Plaque", "Tongue", "Throat"),
                              marginal="nb",
                              zp_cutoff=1,
                              ncores=2)

NB_sim <- simulate_count_scDesign2(model_params=NB_fit,
                                        n_cell_new = 200,
                                        sim_method="copula",
                                        cell_type_prop = cell_type_prop)

saveRDS(NB_sim, file="experiment/scDesign/HMP/NB_sim.rds")


ZINB_fit <- fit_model_scDesign2(data_mat = counts_mat,
                                cell_type_sel=c("Plaque", "Tongue", "Throat"),
                                marginal="zinb",
                                zp_cutoff=1,
                                ncores=2)


ZINB_sim <- simulate_count_scDesign2(model_params=ZINB_fit,
                                   n_cell_new = 200,
                                   sim_method="copula",
                                   cell_type_prop = cell_type_prop)

saveRDS(ZINB_sim, file="experiment/scDesign/HMP/ZINB_sim.rds")


full_fit <- fit_model_scDesign2(data_mat = counts_mat,
                                cell_type_sel=c("Plaque", "Tongue", "Throat"),
                                marginal="auto_choose",
                                zp_cutoff=1,
                                ncores=2)

full_sim <- simulate_count_scDesign2(model_params=ZINB_fit,
                                     n_cell_new = 200,
                                     sim_method="copula",
                                     cell_type_prop = cell_type_prop)

saveRDS(full_sim, file="experiment/scDesign/HMP/full_sim.rds")


saveRDS(list(poisson=poisson_fit, NB=NB_fit,
             ZINB=ZINB_fit, full=full_fit),
        file="experiment/scDesign/HMP/fitted_parameters.rds")



