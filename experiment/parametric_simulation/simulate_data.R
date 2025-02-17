
rm(list=ls())
library(SpiecEasi)
# marginal distribution parameters
taxa_distributions <- read.csv("experiment/parametric_simulation/taxa_parameters.csv",
                               row.names=1)

nsample <- 200
ntaxa <- 300

subset_taxa_distributions <- taxa_distributions[sample(nrow(taxa_distributions), ntaxa), ]

# set up correlation structure
graph <- make_graph('cluster', ntaxa, ntaxa)
prec <- graph2prec(graph, targetCondition=10)
Cov <- prec2cov(prec)
Cor   <- cov2cor(Cov)


# simulate counts from poisson distribution
poi_counts <- rmvpois(n=nsample, mu=subset_taxa_distributions$poi_mean,
                      Sigma=Cor)

nb_counts <- rmvnegbin(n=nsample, mu=subset_taxa_distributions$nb_mean,
                      Sigma=Cor, ks=subset_taxa_distributions$nb_size)

zpois_counts <- rmvzipois(n=nsample, lambdas=subset_taxa_distributions$zpois_mean,
                       Sigma=Cor, ps=subset_taxa_distributions$zpois_zprob)

zinb_counts <- rmvzinegbin(n=nsample, ps=subset_taxa_distributions$zinb_zprob,
                           munbs=subset_taxa_distributions$zinb_mean,
                           ks=subset_taxa_distributions$zinb_size,
                           Sigma=Cor)


# TODO: change distribution parameters for setting up DA taxa and different sample types


