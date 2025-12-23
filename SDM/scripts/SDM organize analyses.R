

# Group species by number of observations to define SDMs algorithms
# Pironon et al. (2024). The global distribution of plants used by humans. Science, 383(6680), 293-297.


# 
library(tidyverse)
library(terra)
library(rnaturalearth)

# 
# setwd('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity')


# observations
load("data/data_final.RData")
scl_occurrences <- data_final$occurrences



# Divide species in three groups: raw (max 2 obs), IDW (3-9 obs), ensemble (min 10 obs)
scl_sdm_alg <- table(scl_occurrences$scientific_name) %>% as.data.frame()
colnames(scl_sdm_alg) <- c('species', 'n_obs')

# order species by occurrences
scl_sdm_alg <- scl_sdm_alg[order(scl_sdm_alg$n_obs, decreasing=F),]

# add algorithm to run
scl_sdm_alg$algorithm <- NA
scl_sdm_alg$algorithm[scl_sdm_alg$n_obs<=2] <- 'raw'
scl_sdm_alg$algorithm[scl_sdm_alg$n_obs>2 & scl_sdm_alg$n_obs<=14] <- 'IDW'
scl_sdm_alg$algorithm[scl_sdm_alg$n_obs>=15 & scl_sdm_alg$n_obs<=49] <- 'maxent'
scl_sdm_alg$algorithm[scl_sdm_alg$n_obs>=50] <- 'ensemble'

table(scl_sdm_alg$algorithm)


# EOO for paper
scl_sdm_alg$EOO <- NA
for (i in 1:nrow(scl_sdm_alg)) {
  scl_sdm_alg$EOO[i] <- data_final$occurrences[data_final$occurrences$scientific_name==scl_sdm_alg$species[i],c('x','y')] %>%
    vect(geom=c('x','y'), crs='epsg:4326') %>% convHull() %>% expanse(unit='km')
  print(i)
}
scl_sdm_alg %>% group_by(algorithm) %>% summarise(mean=mean(EOO), sd=sd(EOO))


# simplified species name
scl_sdm_alg$species_simp <- sub(' ', '_', scl_sdm_alg$species %>% stringr::word(1, 2))


# Set modelling area (AOI): ecoregions with any Scleria observations

# import predictors
# unc_rst_predictors_2.5min <- terra::rast('SDM/unc_rst_predictors_2.5min.tiff')

# import ecoregions
# AOI <- ecoregions %>% terra::intersect(vect(data_final[['occurrences']], geom=c('x','y'), 'epsg:4326'))
# AOI <- AOI$ECO_NAME %>% unique()
# AOI <- terra::subset(ecoregions, ecoregions$ECO_NAME %in% AOI)
# writeVector(AOI, 'AOI/AOI.shp')
AOI <- terra::vect('SDM/AOI/AOI.shp')

# crop predictors and save
# crop_unc_rst_predictors_2.5min <- unc_rst_predictors_2.5min %>% terra::crop(terra::aggregate(AOI), mask=T)
# writeRaster(crop_unc_rst_predictors_2.5min, 'SDM/crop_unc_rst_predictors_2.5min.tiff')
crop_unc_rst_predictors_2.5min <- terra::rast('SDM/crop_unc_rst_predictors_2.5min.tiff')



# coastlines and countries for plots
coastlines <- ne_coastline(scale = 110, returnclass = "sf") %>% vect()
countries <- ne_countries(scale = 110, type = "countries", returnclass = "sf") %>% vect()


