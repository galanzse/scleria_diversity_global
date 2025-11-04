

# CALCULATE MORAN'S I PER SPECCIES USING PRESENCE DATA FOR PAPER
library(tidyverse)
library(terra)
library(ape) # moran.i
library(geosphere) # distm


# check data format
occurrences <- pts_observed %>% terra::project('epsg:4326') %>% geom()
occurrences <- occurrences[,c('x','y')] %>% as.data.frame()
occurrences$species <- pts_observed$epithet

# project raster to equal area
raster <- terra::rast('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/worldclim/wc2.1_30s historical_19702000/wc2.1_30s_bio_1.tif')

# loop
IMoran_scleria <- vector(length=length(unique(occurrences$species)))
names(IMoran_scleria) <- unique(occurrences$species)

for (r in 1:length(IMoran_scleria)) {
    
  # calculate Moran's I
  temp <- occurrences[occurrences$species==names(IMoran_scleria[r]),]
  colnames(temp)[colnames(temp)=='lon'] <- 'x'; colnames(temp)[colnames(temp)=='lat'] <- 'y'
  temp <- temp[,c('x','y')]

  # handle sample size
  if (nrow(temp)<3) next
  
  # extract values from the predictor variable for the points of presence
  xy.variables <- na.omit(
    data.frame(temp,
               terra::extract(x=raster, y=temp, df=TRUE, cells=FALSE, ID=FALSE)
    )
  )
  
  # rename variable to call it later
  colnames(xy.variables)[!(colnames(xy.variables)%in%c('x','y'))] <- 'var'
  
  # matrix of inverse distances between occurrences
  xy.distancias <- 1/ (geosphere::distm(x=xy.variables[, c("x","y")], fun=distGeo) / 1000)
  xy.distancias[!is.finite(xy.distancias)] <- 0 # replace Inf with zero
  diag(xy.distancias) <- 0 # diagonal to 0
  
  # Moran's I 
  IMoran_scleria[r] <- ape::Moran.I(xy.variables$var, xy.distancias, na.rm=TRUE)$observed

  # progress
  print(paste('-- ', round(r/length(IMoran_scleria),2)*100, '% --' , sep=''))
  
}

# save
IMoran_scleria[IMoran_scleria==0] <- NA
mean(IMoran_scleria, na.rm=T) # 0.25
sd(IMoran_scleria, na.rm=T) # 0.39


