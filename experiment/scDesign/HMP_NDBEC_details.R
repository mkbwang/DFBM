
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


set_truth <- function(distribution_truth, metadata, count_mat){
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

  return(normalized_true_abundance)

}

parameter_truth <- readRDS("experiment/scDesign/HMP/fitted_parameters.rds")
zinb_truth <- parameter_truth$ZINB
metadata_zinb <- read.csv("experiment/scDesign/HMP/metadata_zinb.csv")
count_mat_zinb <- read.table("experiment/scDesign/HMP/zinb_sim.txt",
                             header=TRUE, row.names=1, sep="\t") |> as.matrix()


zinb_true_abundance <- set_truth(zinb_truth, metadata_zinb, count_mat_zinb)

ndbec_performance <- function(normalized_true_abundance, ndbec_abundance){
  normalized_ndbec <- row_normalize(ndbec_abundance, libsize=1000)
  performance <- rep(0, nrow(normalized_true_abundance))
  for (j in 1:nrow(normalized_true_abundance)){
    performance[j] <- bray_curtis_distance(normalized_true_abundance[j,],
                                                 normalized_ndbec[j,])
  }
  return(performance)
}

zinb_1 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_1.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_1 <- ndbec_performance(zinb_true_abundance, zinb_1)

zinb_2 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_2.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_2 <- ndbec_performance(zinb_true_abundance, zinb_2)

zinb_3 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_3.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_3 <- ndbec_performance(zinb_true_abundance, zinb_3)

zinb_4 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_4.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_4 <- ndbec_performance(zinb_true_abundance, zinb_4)

zinb_5 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_5.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_5 <- ndbec_performance(zinb_true_abundance, zinb_5)

zinb_6 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_6.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_6 <- ndbec_performance(zinb_true_abundance, zinb_6)


zinb_7 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_7.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_7 <- ndbec_performance(zinb_true_abundance, zinb_7)


zinb_8 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_8.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_8 <- ndbec_performance(zinb_true_abundance, zinb_8)


zinb_9 <- read.table("experiment/scDesign/HMP/NDBEC/zinb_denoise_9.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_zinb_9 <- ndbec_performance(zinb_true_abundance, zinb_9)

combined_performance_zinb <- data.frame(BCDist=c(performance_zinb_1 , performance_zinb_2,
                                                 performance_zinb_3 , performance_zinb_4,
                                                 performance_zinb_5 , performance_zinb_6,
                                                 performance_zinb_7),
                                        Thresholds=rep(c(3,4,5,6,8,11,14), each=200))
combined_performance_zinb$Distribution <- "Zero Inflated Negative Binomial"


nb_truth <- parameter_truth$NB
metadata_nb <- read.csv("experiment/scDesign/HMP/metadata_nb.csv")
count_mat_nb <- read.table("experiment/scDesign/HMP/NB_sim.txt",
                             header=TRUE, row.names=1, sep="\t") |> as.matrix()

nb_true_abundance <- set_truth(nb_truth, metadata_nb, count_mat_nb)


nb_1 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_1.txt",
                   header=TRUE, row.names=1, sep="\t")
performance_nb_1 <- ndbec_performance(nb_true_abundance, nb_1)


nb_2 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_2.txt",
                   header=TRUE, row.names=1, sep="\t")
performance_nb_2 <- ndbec_performance(nb_true_abundance, nb_2)


nb_3 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_3.txt",
                   header=TRUE, row.names=1, sep="\t")
performance_nb_3 <- ndbec_performance(nb_true_abundance, nb_3)

nb_4 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_4.txt",
                   header=TRUE, row.names=1, sep="\t")
performance_nb_4 <- ndbec_performance(nb_true_abundance, nb_4)

nb_5 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_5.txt",
                   header=TRUE, row.names=1, sep="\t")
performance_nb_5 <- ndbec_performance(nb_true_abundance, nb_5)

nb_6 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_6.txt",
                   header=TRUE, row.names=1, sep="\t")
performance_nb_6 <- ndbec_performance(nb_true_abundance, nb_6)

nb_7 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_7.txt",
                     header=TRUE, row.names=1, sep="\t")
performance_nb_7 <- ndbec_performance(nb_true_abundance, nb_7)

# nb_8 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_8.txt",
#                    header=TRUE, row.names=1, sep="\t")
# performance_nb_8 <- ndbec_performance(nb_true_abundance, nb_8)
#
#
# nb_9 <- read.table("experiment/scDesign/HMP/NDBEC/nb_denoise_9.txt",
#                    header=TRUE, row.names=1, sep="\t")
# performance_nb_9 <- ndbec_performance(nb_true_abundance, nb_9)

combined_performance_nb <- data.frame(BCDist=c(performance_nb_1,
                                               performance_nb_2,
                                               performance_nb_3,
                                               performance_nb_4,
                                               performance_nb_5,
                                               performance_nb_6,
                                               performance_nb_7),
                                      Thresholds=rep(c(3,4,5,6,8,11,14), each=200))
combined_performance_nb$Distribution <- "Negative Binomial"


combined_performance <- rbind(combined_performance_nb,
                              combined_performance_zinb)

combined_performance$Thresholds <- factor(combined_performance$Thresholds)

library(ggplot2)
performance_plot <- ggplot(combined_performance, aes(x=Thresholds, y=BCDist)) +
  geom_boxplot()+ labs(x="Number of Thresholds", y="Bray-Curtis Distances")+
  scale_y_continuous(breaks = seq(0.1, 0.4, by = 0.05), limits=c(0.1, 0.4)) +
  facet_grid(cols=vars(Distribution))
