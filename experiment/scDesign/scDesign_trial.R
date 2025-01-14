library(scDesign2)
library(SingleCellExperiment)
library(ggplot2)
library(DuoClustering2018)
library(scran)
# library(tidyverse)


# 10x genomics
sce <- get("sce_filteredExpr10_Zhengmix4eq")(metadata = FALSE)
cell_metadata <- data.frame(colData(sce))
phenotypes <- cell_metadata$phenoid

# colData(sce)$cell_type = as.factor(colData(sce)$phenoid)

# find 500 most variable genes
ngene <- 500
count_mat <- counts(sce)
logcounts(sce) <- log1p(counts(sce))
temp_sce <- modelGeneVar(sce)
chosen <- getTopHVGs(temp_sce, n = ngene)
sce_subset <- sce[chosen,]
counts_subset <- counts(sce_subset)
colnames(counts_subset) <- phenotypes

RNGkind("L'Ecuyer-CMRG")



copula_result <- fit_model_scDesign2(counts_subset,
                                     'b.cells',
                                     sim_method = 'copula',
                                     marginal='auto_choose',
                                     zp_cutoff=1)

set.seed(2024)
full_copula_result <- fit_model_scDesign2(counts_subset,
                                          cell_type_sel = c("b.cells", "cd14.monocytes", "naive.cytotoxic", "regulatory.t"),
                                          sim_method = 'copula',
                                          marginal='auto_choose',
                                          zp_cutoff=1, ncores=4)

saveRDS(full_copula_result,
        file="~/UM/Research/MDAWG/NDBEC/experiment/scDesign/10xgenomics_2017_params.rds")

cell_proportions <- table(phenotypes) / length(phenotypes)
simulated_counts <- simulate_count_scDesign2(model_params=full_copula_result,
                                             n_cell_new = length(phenotypes),
                                             cell_type_prop=cell_proportions)

write.table(data.frame(t(simulated_counts)), sep='\t',
          file = "~/UM/Research/MDAWG/NDBEC/experiment/scDesign/10xgenomics_simulated_count.tsv",
          quote=F)

# set.seed(123)
# example_simu <- scdesign3(
#   sce = sce,
#   assay_use = "counts",
#   celltype = "cell_type",
#   pseudotime = NULL,
#   spatial = NULL,
#   other_covariates = NULL,
#   mu_formula = "cell_type",
#   sigma_formula = "cell_type",
#   family_use = "nb",
#   n_cores = 2,
#   usebam = FALSE,
#   corr_formula = "cell_type",
#   copula = "gaussian",
#   DT = TRUE,
#   pseudo_obs = FALSE,
#   return_model = FALSE,
#   nonzerovar = FALSE,
#   parallelization = "pbmcmapply"
# )


