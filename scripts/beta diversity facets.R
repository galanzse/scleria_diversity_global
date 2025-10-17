

# CALCULATE BETA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY


library(tidyverse)
library(terra)
library(vegan)
library(BAT)
library(ade4)


# load data for analyses
source('scripts/import and prepare data.R')



# Presence matrix
expected_assemblages[,-c(1:3)]



# Taxonomic beta diversity: Bray-Curtis
scl_beta_tax <- vegdist(expected_assemblages[,-c(1:3)], method="bray", binary=FALSE, diag=FALSE, upper=TRUE)

# Functional beta diversity: Btotal (total) = Brepl (replacement) + Brich (species loss/gain)
scl_beta_fun <- BAT::beta(expected_assemblages[,-c(1:3)], scl_dendrogram, func="jaccard", abund=F, raref=0, comp=F)

# Phylogenetic beta diversity
scl_beta_phy <- BAT::beta(expected_assemblages[,-c(1:3)], scl_phylogeny, func="jaccard", abund=F, raref=0, comp=F)

# save results
scl_beta_div <- list()
scl_beta_div[['taxonomic']] <- scl_beta_tax
scl_beta_div[['functional']] <- scl_beta_fun
scl_beta_div[['phylogenetic']] <- scl_beta_phy

save(scl_beta_div, file='results/scl_beta_div.Rdata')



# correlation between matrices
mantel.rtest


