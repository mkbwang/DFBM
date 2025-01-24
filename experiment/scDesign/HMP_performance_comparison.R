
rm(list=ls())
library(vegan)

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
  data <- rbind(u, v)
  distance <- vegdist(data, method = "bray")
  return(distance)
}


performance_summary <- function(distribution_truth, metadata, count_mat,
                                ndbec_denoise, magic_denoise, saver_denoise, scimpute_denoise,
                                dca_denoise, autoclass_denoise, alra_denoise){

  type_counts <- table(metadata$Type) # counts of different sample types

  # calculate the expected abundance of each taxa for different sample types
  plaque_truth <- distribution_truth$Plaque
  plaque_mean <- rep(0, ncol(count_mat))
  plaque_mean[setdiff(1:ncol(count_mat), plaque_truth$gene_sel3)] <-
    plaque_truth$marginal_param1[, 3] * (1-plaque_truth$marginal_param1[, 1])

  tongue_truth <- distribution_truth$Tongue
  tongue_mean <- rep(0, ncol(count_mat))
  tongue_mean[setdiff(1:ncol(count_mat), tongue_truth$gene_sel3)] <-
    tongue_truth$marginal_param1[, 3] * (1-tongue_truth$marginal_param1[, 1])


  throat_truth <- distribution_truth$Throat
  throat_mean <- rep(0, ncol(count_mat))
  throat_mean[setdiff(1:ncol(count_mat), throat_truth$gene_sel3)] <-
    throat_truth$marginal_param1[, 3] * (1-throat_truth$marginal_param1[, 1])

  ## calculate true normalized abundance
  truth_abundance_mat <- rbind(t(replicate(type_counts['Plaque'], plaque_mean)),
                               t(replicate(type_counts['Tongue'], tongue_mean)),
                               t(replicate(type_counts['Throat'], throat_mean)))
  normalized_true_abundance <- row_normalize(truth_abundance_mat, libsize=1000)


  ## calculate the normalized count matrix
  normalized_count_mat <- row_normalize(count_mat, libsize=1000)

  ## bray curtis distance between observed counts and true abundances
  raw_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    raw_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                               normalized_count_mat[j,])
  }
  raw_performance_df <- data.frame(BCdist=raw_performance, Source="Raw",
                                   Type=metadata$Type)


  ##  performance of NDBEC
  normalized_ndbec <- row_normalize(ndbec_denoise, libsize=1000)

  ndbec_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    ndbec_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                                 normalized_ndbec[j,])
  }
  ndbec_performance_df <- data.frame(BCdist=ndbec_performance, Source="NDBEC",
                                     Type=metadata$Type)


  ## performance of MAGIC
  magic_denoise_sq <- magic_denoise^2
  normalized_magic <- row_normalize(magic_denoise_sq, libsize=1000)
  magic_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    magic_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                                 normalized_magic[j,])
  }
  magic_performance_df <- data.frame(BCdist=magic_performance, Source="MAGIC",
                                     Type=metadata$Type)


  ## performance of SAVER
  normalized_saver <- row_normalize(saver_denoise, libsize=1000)
  saver_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    saver_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                                 normalized_saver[j,])
  }
  saver_performance_df <- data.frame(BCdist=saver_performance, Source="SAVER",
                                     Type=metadata$Type)


  ## performance of scImpute
  normalized_scimpute <- row_normalize(t(scimpute_denoise), libsize=1000)
  scimpute_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    scimpute_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                                 normalized_scimpute[j,])
  }
  scimpute_performance_df <- data.frame(BCdist=scimpute_performance, Source="scImpute",
                                     Type=metadata$Type)

  ## performance of ALRA
  normalized_alra <- row_normalize(alra_denoise, libsize=1000)
  alra_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    alra_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                                normalized_alra[j,])
  }
  alra_performance_df <- data.frame(BCdist=alra_performance, Source="ALRA",
                                    Type=metadata$Type)


  ## performance of DCA
  dca_denoise <- t(dca_denoise)
  normalized_dca <- row_normalize(dca_denoise, libsize=1000)
  dca_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    dca_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                               normalized_dca[j,])
  }
  dca_performance_df <- data.frame(BCdist=dca_performance, Source="DCA",
                                   Type=metadata$Type)

  ## performance of AutoClass
  normalized_autoclass <- row_normalize(autoclass_denoise, libsize=1000)
  autoclass_performance <- rep(0, nrow(count_mat))
  for (j in 1:nrow(count_mat)){
    autoclass_performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                               normalized_autoclass[j,])
  }
  autoclass_performance_df <- data.frame(BCdist=autoclass_performance, Source="AutoClass",
                                   Type=metadata$Type)


  combined_performance <- rbind(raw_performance_df,
                                magic_performance_df,
                                saver_performance_df,
                                scimpute_performance_df,
                                dca_performance_df,
                                autoclass_performance_df,
                                alra_performance_df,
                                ndbec_performance_df)

  combined_performance$Source <- factor(combined_performance$Source,
                                        levels=c("Raw", "MAGIC", "SAVER", "scImpute", "DCA",
                                                 "AutoClass", "ALRA", "NDBEC"))

  return(combined_performance)

}

library(ggplot2)
plot_summary <- function(summary_df, min=0, max=0.8, gap=0.1){

  summary_df$cm <- 1
  summary_df$cm[summary_df$Source == "Raw"] <- 2
  summary_df$cm[summary_df$Source == "NDBEC"] <- 3
  summary_df$cm <- factor(summary_df$cm)

  cmaps <- c("#575656", "black", "#105aa3")

  output_plot <- ggplot(summary_df, aes(x=Source, y=BCdist, color=cm)) +
    geom_boxplot()+ labs(x=NULL, y="Bray-Curtis Distances") +
    scale_color_manual(values=cmaps)+
    scale_y_continuous(breaks = seq(min, max, by = gap), limits=c(min, max)) +
    theme(legend.position = "none")+facet_grid(cols=vars(Distribution))

  return(output_plot)

}


parameter_truth <- readRDS("experiment/scDesign/HMP/fitted_parameters.rds")


# ZINB
zinb_truth <- parameter_truth$ZINB
count_mat_zinb <- read.table("experiment/scDesign/HMP/zinb_sim.txt",
                        header=TRUE, row.names=1, sep="\t") |> as.matrix()
metadata_zinb <- read.csv("experiment/scDesign/HMP/metadata_zinb.csv")
ndbec_zinb <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_4.txt",
                            header=TRUE, row.names=1, sep="\t")
magic_zinb <- read.table("experiment/scDesign/HMP/MAGIC/magic_zinb.tsv",
                         header=TRUE, row.names=1, sep="\t")
saver_zinb <- read.table("experiment/scDesign/HMP/SAVER/SAVER_zinb.tsv",
                            header=TRUE, row.names=1, sep="\t")
scimpute_zinb <- read.table("experiment/scDesign/HMP/scImpute/scImpute_zinb.tsv",
                            header=TRUE, row.names=1, sep=" ")
autoclass_zinb <- read.table("experiment/scDesign/HMP/AutoClass/AutoClass_zinb.tsv",
                             header=TRUE, row.names=1, sep="\t")
alra_zinb <- read.table("experiment/scDesign/HMP/ALRA/ALRA_zinb.tsv",
                           header=TRUE, row.names=1, sep="\t")
dca_zinb <- read.table("experiment/scDesign/HMP/DCA/DCA_zinb.tsv",
                          header=TRUE, row.names=1, sep="\t")

zinb_summary <- performance_summary(zinb_truth, metadata_zinb, count_mat_zinb,
                                    ndbec_zinb, magic_zinb, saver_zinb, scimpute_zinb,
                                    dca_zinb, autoclass_zinb, alra_zinb)
zinb_summary$Distribution <- "Zero Inflated Negative Binomial Distribution"



plot_zinb <- plot_summary(zinb_summary, min=0, max=0.8, gap=0.1)



# NB
nb_truth <- parameter_truth$NB
count_mat_nb <- read.table("experiment/scDesign/HMP/NB_sim.txt",
                           header=TRUE, row.names=1, sep="\t") |> as.matrix()
metadata_nb <- read.csv("experiment/scDesign/HMP/metadata_nb.csv")
ndbec_nb <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_4.txt",
                         header=TRUE, row.names=1, sep="\t")
magic_nb <- read.table("experiment/scDesign/HMP/MAGIC/magic_nb.tsv",
                         header=TRUE, row.names=1, sep="\t")
saver_nb <- read.table("experiment/scDesign/HMP/SAVER/SAVER_nb.tsv",
                         header=TRUE, row.names=1, sep="\t")
scimpute_nb <- read.table("experiment/scDesign/HMP/scImpute/scImpute_nb.tsv",
                            header=TRUE, row.names=1, sep=" ")
alra_nb <- read.table("experiment/scDesign/HMP/ALRA/ALRA_nb.tsv",
                        header=TRUE, row.names=1, sep="\t")
dca_nb <- read.table("experiment/scDesign/HMP/DCA/DCA_nb.tsv",
                       header=TRUE, row.names=1, sep="\t")
autoclass_nb <- read.table("experiment/scDesign/HMP/AutoClass/AutoClass_nb.tsv",
                             header=TRUE, row.names=1, sep="\t")
nb_summary <- performance_summary(nb_truth, metadata_nb, count_mat_nb,
                                    ndbec_nb, magic_nb, saver_nb, scimpute_nb,
                                    dca_nb, autoclass_nb, alra_nb)
nb_summary$Distribution <- "Negative Binomial Distribution"


combined_summary <- rbind(zinb_summary, nb_summary)

cmaps <- c("#575656", "black", "#105aa3")
overall_performance <- plot_summary(combined_summary, min=0, max=0.8, gap=0.1)
plaque_performance <- plot_summary(combined_summary[combined_summary$Type == "Plaque",],
                                  min=0, max=0.8, gap=0.1)

tongue_performance <- plot_summary(combined_summary[combined_summary$Type == "Tongue",],
                                  min=0, max=0.8, gap=0.1)


throat_performance <- plot_summary(combined_summary[combined_summary$Type == "Throat",],
                                   min=0, max=0.8, gap=0.1)

plot_nb <- plot_summary(nb_summary, min=0, max=0.8, gap=0.1)



# Poisson
# poisson_truth <- parameter_truth$poisson
# count_mat_poisson <- read.table("experiment/scDesign/HMP/poisson_sim.txt",
#                            header=TRUE, row.names=1, sep="\t") |> as.matrix()
# metadata_poisson <- read.csv("experiment/scDesign/HMP/metadata_poisson.csv")
# ndbec_poisson <- read.table("experiment/scDesign/HMP/NDBEC/poisson_denoise_9.txt",
#                        header=TRUE, row.names=1, sep="\t")
# magic_poisson <- read.table("experiment/scDesign/HMP/MAGIC/magic_poisson.tsv",
#                        header=TRUE, row.names=1, sep="\t")
# saver_poisson <- read.table("experiment/scDesign/HMP/SAVER/SAVER_poisson.tsv",
#                        header=TRUE, row.names=1, sep="\t")
# alra_poisson <- read.table("experiment/scDesign/HMP/ALRA/ALRA_poisson.tsv",
#                       header=TRUE, row.names=1, sep="\t")
# dca_poisson <- read.table("experiment/scDesign/HMP/DCA/DCA_poisson.tsv",
#                      header=TRUE, row.names=1, sep="\t")
# poisson_summary <- performance_summary(poisson_truth, metadata_poisson, count_mat_poisson,
#                                   ndbec_poisson, magic_poisson, saver_poisson,
#                                   dca_poisson, alra_poisson)
# plot_poisson <- plot_summary(poisson_summary)




# full
# full_truth <- parameter_truth$full
#
#
# library(dplyr)
# combined_performance_summary <- combined_performance %>%
#   group_by(Source) %>%
#   summarize(mean=mean(BCdist), sd=sd(BCdist))
#
# combined_performance_summary$mean <- round(combined_performance_summary$mean, 3)
# combined_performance_summary$sd <- round(combined_performance_summary$sd, 3)
#
#
# write.csv(combined_performance_summary, "experiment/scDesign/combined_performance_summary.csv")
#
# library(ggplot2)
























