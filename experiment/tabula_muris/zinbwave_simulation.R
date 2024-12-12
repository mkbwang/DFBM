
library(Matrix)
rm(list=ls())
count_mat <- readMM("experiment/tabula_muris/ZINBwave/mouse_fat_template_mat.mtx")
gene_names <- readLines("experiment/tabula_muris/ZINBwave/genes_template.txt")
cell_names <- readLines("experiment/tabula_muris/ZINBwave/cells_template.txt")
rownames(count_mat) <- gene_names
colnames(count_mat) <- cell_names

count_mat_df <- as.data.frame(as.matrix(t(count_mat))) 
write.table(count_mat_df, "experiment/tabula_muris/ZINBwave/template_raw.tsv",
            sep='\t', quote=F)


metadata <- read.table("experiment/tabula_muris/ZINBwave/cells_metadata.tsv",
                       header=T, sep='\t')


metadata$Bcell <- 1*(metadata$cell_ontology_class == "B cell")
metadata$Tcell <- 1*(metadata$cell_ontology_class == "T cell")
metadata$endothelial <- 1*(metadata$cell_ontology_class == "endothelial cell")
metadata$mesenchymal <- 1*(metadata$cell_ontology_class == "mesenchymal stem cell of adipose")
metadata$myeloid <- 1*(metadata$cell_ontology_class == "myeloid cell")

metadata_subset <- metadata[, c("Bcell", "Tcell", "endothelial", "mesenchymal", "myeloid")]

# visualize the original expressions on the log scale
log_count_mat <- log10(as.matrix(count_mat) + 1)
hist(sqrt(as.matrix(count_mat)), nclass=20)
rownames(log_count_mat) <- NULL
colnames(log_count_mat) <- NULL
source("experiment/visualization_utils.R")
viz_original <- plot_heatmap(log_count_mat, legend_title="Log10Count")

# Set up singlecellexperiment object
library(SingleCellExperiment)
library(zinbwave)


sce <- SingleCellExperiment(list(counts=as.matrix(count_mat)),
                            colData = metadata_subset)  

library(BiocParallel)
zinb_fit <- zinbFit(Y = sce, X = "~Bcell+Tcell+endothelial+mesenchymal+myeloid-1",
                    K=0, maxiter.optimize=40, which_assay="counts",
                    BPPARAM=MulticoreParam(4))

X_mat <- zinb_fit@X
beta_hat_mu <- zinb_fit@beta_mu
beta_hat_pi <- zinb_fit@beta_pi
V_mat <- t(zinb_fit@V)
gamma_hat_mu <- t(zinb_fit@gamma_mu)
gamma_hat_pi <- t(zinb_fit@gamma_pi)

log_mu_hat <- X_mat %*% beta_hat_mu + gamma_hat_mu %*% V_mat
mu_hat <- exp(log_mu_hat) 

viz_logmu <- plot_heatmap(log10(t(mu_hat)+1), legend_title="Log10mu")


mu_hat_df <- data.frame(mu_hat)
rownames(mu_hat_df) <- colnames(count_mat)
colnames(mu_hat_df) <- rownames(count_mat)
write.table(mu_hat_df, "experiment/tabula_muris/ZINBwave/mu_hat.tsv",
            sep='\t', quote=F)


logit_pi_hat <- X_mat %*% beta_hat_pi + gamma_hat_pi %*% V_mat
pi_hat <- exp(logit_pi_hat)/(1+exp(logit_pi_hat))
viz_pi <- plot_heatmap(t(1-pi_hat), legend_title="NonzeroProb",
                       min=0, max=1, colormap="inferno")


zinb_mu_hat <- mu_hat * (1-pi_hat)
log_zinb_mu_hat <- log10(zinb_mu_hat+1)
viz_zinb_mu <- plot_heatmap(t(log_zinb_mu_hat), legend_title="log10ZINB",
                            min=0, max=3)

pi_hat_df <- data.frame(pi_hat)
rownames(pi_hat_df) <- colnames(count_mat)
colnames(pi_hat_df) <- rownames(count_mat)
write.table(pi_hat_df, "experiment/tabula_muris/ZINBwave/pi_hat_sample.tsv",
            sep='\t', quote=F)


theta_hat <- exp(zinb_fit@zeta)
theta_hat_df <- data.frame(theta_hat)
rownames(theta_hat_df) <- rownames(count_mat)
write.table(theta_hat_df, "experiment/tabula_muris/ZINBwave/theta_hat.tsv",
            sep='\t', quote=F)


## simulate count data based on the estimate pi, mu and theta

new_counts <- matrix(0, nrow=ncol(count_mat), 
                     ncol=nrow(count_mat))


for (i in 1:ncol(count_mat)) {
  for (j in 1:nrow(count_mat)) {
    if (runif(1) > pi_hat[i, j]) {
      new_counts[i, j] <- rnbinom(1, mu = mu_hat[i, j], size = theta_hat[j])
    }
  }
}
log_new_counts <- log10(new_counts + 1)


viz_logcounts <- plot_heatmap(t(log_new_counts),
                              legend_title="Log10Count")

new_counts_df <- data.frame(new_counts)
rownames(new_counts_df) <- colnames(count_mat)  
colnames(new_counts_df) <- rownames(count_mat)  

write.table(new_counts_df, "experiment/tabula_muris/ZINBwave/simulated_counts.tsv",
            sep='\t', quote=F)



