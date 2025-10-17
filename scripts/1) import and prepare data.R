

# PREPARE DATA FOR THE CALCULATION OF ALPHA AND BETA DIVERSITY, INCLUDING NULL MODELS


library(tidyverse)
library(rnaturalearth)
library(terra)
library(cluster) # daisy
library(picante) # randomizeMatrix



# data
load("data/data_final.RData")

# infrageneric classification, treat sections as genera for trait imputation
scl_taxonomy <- data_final[['taxa']] %>% filter(!(is.na(section)))
# simplify names to control for automatic changes
scl_taxonomy$epithet <- sapply(strsplit(scl_taxonomy$scientific_name, " "), function(x) x[2])
# change (-) to avoid posterior errors in table reading
scl_taxonomy$epithet[scl_taxonomy$epithet=='baroni-clarkei'] <- 'baroniclarkei'
scl_taxonomy$epithet[scl_taxonomy$epithet=='flagellum-nigrorum'] <- 'flagellumnigrorum'
scl_taxonomy$epithet[scl_taxonomy$epithet=='novae-hollandiae'] <- 'novaehollandiae'



# OCCURRENCES ####
scl_occurrences <- data_final[['occurrences']]
scl_occurrences <- scl_occurrences %>%
  merge(scl_taxonomy[,c('scientific_name','epithet')], by='scientific_name') %>%
  dplyr::select(epithet, x, y)


# world_lines
world_lines <- ne_countries(scale=10, type="countries", continent=NULL, country=NULL, geounit=NULL, sovereignty=NULL, returnclass="sf") %>% terra::vect() %>% terra::project('+proj=eqearth')

# world_grid
world_grid <- terra::rast('results/maps/scleria_richness_50km.tiff')
world_grid[!is.na(world_grid)] <- 0
# world_grid %>% cellSize(unit='km') %>% sqrt()



# TRAITS ####
source('scripts/1.1) traits.R')



# PHYLOGENY ####
load("results/imputed_trees.RData")

# # plot
# scl_tree1 <- imputed_trees[[15]]
# 
# # Pick colors for each group
# group_colors <- c('Scleria'='forestgreen', 'Hypoporum'='orange', 'Trachylomia'='violet', 'Browniae'='blue')
# 
# # Map species to colors
# tip_colors <- group_colors[scl_taxonomy$subgenus]
# 
# # Make sure the order matches tree$tip.label
# names(tip_colors) <- scl_taxonomy$epithet
# tip_colors <- tip_colors[scl_tree1$tip.label]
# plot.phylo(scl_tree1, tip.color=tip_colors, cex=0.6, main='Phylogeny')
# legend("bottomleft", y.intersp=0.3, x.intersp=0.5,
#        legend=names(group_colors), col=group_colors,
#        pch=19, pt.cex=1, bty="n",  title="")



# DEFINE ASSEMBLAGES ####


# # observed assemblages
# pts_observed <- terra::vect(scl_occurrences, geom=c('x','y'), crs='epsg:4326') %>%
#   terra::project(world_grid)
# 
# observed_assemblages <- world_grid %>% terra::extract(pts_observed, cells=T, xy=T) %>% dplyr::select(cell, x, y)
# observed_assemblages$epithet <- pts_observed$epithet
# 
# # check: remove species without trait and mollecular data
# observed_assemblages <- observed_assemblages %>% subset(epithet %in% v_species)
# 
# # df into presence matrix
# observed_assemblages$presence <- 1
# observed_assemblages <- observed_assemblages %>% unique() %>%
#   pivot_wider(names_from='epithet', values_from='presence') %>% as.data.frame()
# observed_assemblages[is.na(observed_assemblages)] <- 0
# rownames(observed_assemblages) <- observed_assemblages$cell
# 
# write.table(observed_assemblages, 'results/observed_assemblages.txt')
observed_assemblages <- read.table('results/observed_assemblages.txt')


# # expected occurrences
# expected_assemblages <- terra::rast('results/maps/scleria_probocc_50km.tiff') %>%
#   as.data.frame(cells=TRUE, xy=TRUE)
# 
# # fix names to match other datasets
# colnames(expected_assemblages)[4:ncol(expected_assemblages)] <- sub(".*_", "", colnames(expected_assemblages)[4:ncol(expected_assemblages)])
# colnames(expected_assemblages)[colnames(expected_assemblages)=='baroni-clarkei'] <- 'baroniclarkei'
# colnames(expected_assemblages)[colnames(expected_assemblages)=='flagellum-nigrorum'] <- 'flagellumnigrorum'
# colnames(expected_assemblages)[colnames(expected_assemblages)=='novae-hollandiae'] <- 'novaehollandiae'
# 
# # select species common to all datasets
# expected_assemblages <- expected_assemblages[,c('cell','x','y',v_species)]
# 
# # exclude cells with prob. occurrence <1 (retain cells with 100% accumulated probability of species occurrence)
# expected_assemblages <- expected_assemblages[rowSums(expected_assemblages[,-c(1:3)])>=1,] # ~40000 assemblages species out
# hist(rowSums(expected_assemblages[,-c(1:3)]), 40)
# 
# # round to save space
# expected_assemblages[,4:ncol(expected_assemblages)] <- round(expected_assemblages[,4:ncol(expected_assemblages)], 2)
# 
# # save
# write.table(expected_assemblages, 'results/expected_assemblages.txt')
expected_assemblages <- read.table('results/expected_assemblages.txt')


