
rm(list=ls())

row_normalize <- function(mymat, libsize=600){
  t(apply(mymat, 1, function(row) {
    row_sum <- sum(row)
    row * (libsize / row_sum)
  }))
}
bray_curtis_distance <- function(u, v) {
  if (length(u) != length(v)) {
    stop("The vectors must be of the same length.")
  }
  sum_min <- sum(pmin(u, v))
  sum_total <- sum(u) + sum(v)
  bray_curtis <- 1 - (2 * sum_min / sum_total)
  return(bray_curtis)
}



count_mat <- read.table("experiment/scDesign/10xgenomics_simulated_count.tsv",
                        header=TRUE, row.names=1, sep="\t") |> as.matrix()
normalized_count_mat <- row_normalize(count_mat)

celltypes <- gsub('[0-9]+', '', rownames(count_mat))
celltypes <- gsub("\\.$", "", celltypes)
type_count <- table(celltypes)


parameter_truth <- readRDS("experiment/scDesign/10xgenomics_2017_params.rds")


bcell_truth <- parameter_truth$b.cells
bcells_mean <- rep(0, ncol(count_mat))
bcells_mean[setdiff(1:500, bcell_truth$gene_sel3)] <-
  bcell_truth$marginal_param1[, 3] * (1-bcell_truth$marginal_param1[, 1])





monocyte_truth <- parameter_truth$cd14.monocytes
monocyte_mean <- rep(0, ncol(count_mat))
monocyte_mean[setdiff(1:500, monocyte_truth$gene_sel3)] <-
  monocyte_truth$marginal_param1[,3] * (1-monocyte_truth$marginal_param1[,1])


ncytotoxic_truth <- parameter_truth$naive.cytotoxic
ncytotoxic_mean <- rep(0, ncol(count_mat))
ncytotoxic_mean[setdiff(1:500, ncytotoxic_truth$gene_sel3)] <-
  ncytotoxic_truth$marginal_param1[, 3] * (1-ncytotoxic_truth$marginal_param1[, 1])


rt_truth <- parameter_truth$regulatory.t
rt_mean <- rep(0, ncol(count_mat))
rt_mean[setdiff(1:500, rt_truth$gene_sel3)] <-
  rt_truth$marginal_param1[, 3] * (1-rt_truth$marginal_param1[, 1])

truth_expression_mat <- rbind(t(replicate(type_count['b.cells'], bcells_mean)),
                  t(replicate(type_count['cd.monocytes'], monocyte_mean)),
                  t(replicate(type_count['naive.cytotoxic'], ncytotoxic_mean)),
                  t(replicate(type_count['regulatory.t'], rt_mean)))
normalized_true_expression <- row_normalize(truth_expression_mat)


raw_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  raw_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                               normalized_count_mat[j,])
}

raw_performance_df <- data.frame(BCdist=raw_performance, Source="Raw")


# NDBEC
ndbec_denoise <- read.table("experiment/scDesign/ndbec_trial_expected_counts.tsv",
                            header=TRUE, row.names=1, sep="\t")
normalized_ndbec <- row_normalize(ndbec_denoise)

ndbec_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  ndbec_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                               normalized_ndbec[j,])
}
ndbec_performance_df <- data.frame(BCdist=ndbec_performance, Source="Proposed")



# ALRA
alra_denoise <- read.table("experiment/scDesign/competitors/denoised_counts_ALRA.tsv",
                           header=TRUE, row.names=1, sep="\t")
normalized_alra <- row_normalize(alra_denoise)
alra_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  alra_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                               normalized_alra[j,])
}
alra_performance_df <- data.frame(BCdist=alra_performance, Source="ALRA")



# AutoClass
autoclass_denoise <- read.table("experiment/scDesign/competitors/denoised_counts_autoclass.tsv",
                                header=TRUE, row.names=1, sep="\t")
normalized_autoclass <- row_normalize(autoclass_denoise)
autoclass_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  autoclass_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                              normalized_autoclass[j,])
}
autoclass_performance_df <- data.frame(BCdist=autoclass_performance, Source="AutoClass")


# DCA
dca_denoise <- read.table("experiment/scDesign/competitors/denoised_counts_DCA.tsv",
                          header=TRUE, row.names=1, sep="\t")
dca_denoise <- t(dca_denoise)
normalized_dca <- row_normalize(dca_denoise)
dca_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  dca_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                                   normalized_dca[j,])
}
dca_performance_df <- data.frame(BCdist=dca_performance, Source="DCA")



# MAGIC
magic_denoise <- read.table("experiment/scDesign/competitors/denoised_counts_MAGIC.tsv",
                            header=TRUE, row.names=1, sep="\t")
normalized_magic <- row_normalize(magic_denoise)
magic_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  magic_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                             normalized_magic[j,])
}
magic_performance_df <- data.frame(BCdist=magic_performance, Source="MAGIC")




# SAVER
saver_denoise <- read.table("experiment/scDesign/competitors/denoised_counts_SAVER.tsv",
                            header=TRUE, row.names=1, sep="\t")
normalized_saver <- row_normalize(saver_denoise)
saver_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  saver_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                               normalized_saver[j,])
}
saver_performance_df <- data.frame(BCdist=saver_performance, Source="SAVER")


# scImpute
scimpute_denoise <- read.table("experiment/scDesign/competitors/denoised_counts_scimpute.tsv",
                               header=TRUE, row.names=1, sep=" ") |> t()
normalized_scimpute <- row_normalize(scimpute_denoise)
scimpute_performance <- rep(0, nrow(count_mat))
for (j in 1:nrow(count_mat)){
  scimpute_performance[j] <- bray_curtis_distance(normalized_true_expression[j,],
                                               normalized_scimpute[j,])
}
scimpute_performance_df <- data.frame(BCdist=scimpute_performance, Source="scImpute")





combined_performance <- rbind(raw_performance_df, scimpute_performance_df,
                              magic_performance_df, alra_performance_df,
                              saver_performance_df, autoclass_performance_df,
                              dca_performance_df, ndbec_performance_df)

combined_performance$Source <- factor(combined_performance$Source,
                                      levels=c("Raw", "scImpute", "MAGIC", "ALRA", "AutoClass",
                                               "DCA", "SAVER", "Proposed"))

library(dplyr)
combined_performance_summary <- combined_performance %>%
  group_by(Source) %>%
  summarize(mean=mean(BCdist), sd=sd(BCdist))

combined_performance_summary$mean <- round(combined_performance_summary$mean, 3)
combined_performance_summary$sd <- round(combined_performance_summary$sd, 3)


write.csv(combined_performance_summary, "experiment/scDesign/combined_performance_summary.csv")

library(ggplot2)

combined_performance$cm <- 1
combined_performance$cm[combined_performance$Source == "Raw"] <- 2
combined_performance$cm[combined_performance$Source == "Proposed"] <- 3
combined_performance$cm <- factor(combined_performance$cm)

cmaps <- c("#575656", "black", "#105aa3")

ggplot(combined_performance, aes(x=Source, y=BCdist, color=cm)) +
  geom_boxplot()+ labs(x=NULL, y="Bray-Curtis Distances") +
  scale_color_manual(values=cmaps)+
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1)) +
  theme(legend.position = "none")




















