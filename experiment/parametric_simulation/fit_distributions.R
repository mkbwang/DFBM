
rm(list=ls())
library(SpiecEasi)
library(phyloseq)

# simulate from real data
data(amgut1.filt)
# data(amgut2.filt.phy)

ntaxa <- ncol(amgut1.filt)
nsample <- nrow(amgut1.filt)


candidate_parameters <- list()

for (k in 1:10){
  amgut_rawcounts <- amgut1.filt[sample(nsample, 200), ]
  depths <- rowSums(amgut_rawcounts)
  amgut_normalized <- t(apply(amgut_rawcounts, 1, norm_to_total))
  amgut_scaled <- round(amgut_normalized * min(depths))


  ## need to break down this one to generate differentially abundant taxa from different clusters
  poi_params <- get_comm_params(amgut_scaled, mar=2, distr="pois")
  nb_params <- get_comm_params(amgut_scaled, mar=2, distr="negbin")
  zpois_params <- get_comm_params(amgut_scaled, mar=2, distr="zipois")
  zinb_params <- get_comm_params(amgut_scaled, mar=2, distr="zinegbin")



  poi_mean <- unlist(lapply(poi_params, function(estim) estim$lambda))
  nb_mean <- unlist(lapply(nb_params, function(estim) estim$mu))
  nb_size <- unlist(lapply(nb_params, function(estim) estim$size))
  zpois_zprob <- unlist(lapply(zpois_params, function(estim) estim$pstr0))
  zpois_mean <- unlist(lapply(zpois_params, function(estim) estim$lambda))
  zinb_zprob <- unlist(lapply(zinb_params, function(estim) estim$pstr0))
  zinb_mean <- unlist(lapply(zinb_params, function(estim) estim$munb))
  zinb_size <- unlist(lapply(zinb_params, function(estim) estim$size))

  fitted_parameters <- cbind(poi_mean, nb_mean, nb_size, zpois_zprob, zpois_mean,
                             zinb_zprob, zinb_mean, zinb_size) |> as.data.frame()

  candidate_parameters[[k]] <- fitted_parameters

}

candidate_params_df <- do.call(rbind, candidate_parameters)

library(dplyr)

candidate_params_df <- candidate_params_df %>% filter(zinb_mean < 100)
write.csv(candidate_params_df, "experiment/parametric_simulation/taxa_parameters.csv")

