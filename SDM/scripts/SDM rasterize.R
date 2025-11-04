

# rasterize species nobs=1,2


source('SDM/scripts/SDM organize analyses.R')


# species
scl_sdm_raw <- scl_sdm_alg[scl_sdm_alg$algorithm=='raw',]

# richness map at original resolution
r <- terra::rast('C:/Users/javie/Desktop/world_rasters/elevation/wc2.1_2.5m_elev.tif')
# cellSize(r, unit='km') %>% sqrt() # 4.62 km at equator

# set the background cells in the raster to 0
r[!is.na(r)] <- 0

# loop
for (i in 1:nrow(scl_sdm_raw)) {

  # select species and create spatvector
  df_p <- scl_occurrences %>%
    subset(scientific_name%in%scl_sdm_raw$species[i])
  pts_p <- df_p %>% terra::vect(geom=c('x','y'), crs='+proj=longlat')
  
  # plot in ecoregions that intersect the species EOO
  aoi_eco <- AOI %>% terra::mask(terra::convHull(pts_p)) %>% terra::aggregate(dissolve=TRUE)
  
  # crop r to facilitate visualization
  r2 <- r %>% terra::crop(aoi_eco, mask=T)
  
  # create map
  path1 <- 'C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/maps/'
  png(paste(path1, scl_sdm_raw$species[i], '.png', sep=''), bg="white")
  plot(aoi_eco, main=scl_sdm_raw$species[i])
  points(pts_p, col='green', pch=16, cex=2)
  dev.off()
  
  # rasterize, add to base map and save
  temp_rast <- terra::rasterize(pts_p, r2) %>% terra::merge(r2)
  writeRaster(temp_rast, paste('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/predictions/', scl_sdm_raw$species[i], '.tiff', sep=''), overwrite=TRUE)
  
  # progress
  print(paste(round(i/nrow(scl_sdm_raw), 4)*100, '%'))
}


