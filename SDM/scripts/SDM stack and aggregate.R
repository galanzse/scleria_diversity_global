

# COMPILE RESULTS FOR DOWNSTREAM ANALYSES

# it is not really necessary to project to eqearth because species probabilities are not dependent on cellsize, but facilitates interpretation


source('SDM/scripts/SDM organize analyses.R')


# project to eqearth: original and aggregated ####

# base map
r <- terra::rast('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/elevation/wc2.1_2.5m_elev.tif')
names(r) <- 'base'
r[!is.na(r)] <- 0
cellSize(r, unit='km') %>% sqrt()  

# cells must have same cell size for alpha diversity analyses
r <- r %>% terra::project('+proj=eqearth')
e <- ext(c(-17243958, 17243039, -8393409, 8392926))
r <- terra::crop(r, e)
# cellSize(r, unit='km') %>% sqrt() # 2.73*2.73km (6.76 km2) at equator

# aggregate to 50km for alphadiv analyses
r_agg <- r %>% terra::aggregate(fact=18)
# cellSize(r_agg, unit='km') %>% sqrt() # 49.28*49.28km (2429 km2) at equator


# import SDMs, project, extend, merge and save
temp_files <- list.files('SDM/results/predictions', full.names=TRUE)

# loop
for (i in 1:length(temp_files)) {

  # species name
  s = strsplit(temp_files[i], "/")[[1]][4]
  print(s)

  # merge retains the values of the first SpatRaster in the sequence of arguments
  ee_rast <- temp_files[i] %>% terra::rast() %>% terra::project('+proj=eqearth') %>%
    terra::resample(r, method='bilinear') %>% merge(r)

  # resample at lower resolution for alphadiv analyses
  # mean: weighted average prob. of occurrence within any point in the raster cell
  ee_rast50 <- ee_rast %>% terra::resample(y=r_agg, method='bilinear') %>%
    # avoid the creation of new cells when using mean, 'merge' does not solves the issue
    sum(r_agg, na.rm=FALSE)

  # save
  writeRaster(ee_rast, paste('SDM/results/eqearth_predictions/2.8km_original/', s, sep=''), overwrite=TRUE)
  writeRaster(ee_rast50, paste('SDM/results/eqearth_predictions/50km_agg/', s, sep=''), overwrite=TRUE)

  print(paste(round(i/length(temp_files), 4)*100, '%'))
}



# stack and calculate richness ####


# ORIGINAL

# import extended SDMs
scleria_probocc_2.8km <- list.files('SDM/results/eqearth_predictions/2.8km_original', full.names=TRUE)
# exclude .aux files and prepare stack
scleria_probocc_2.8km <- scleria_probocc_2.8km[seq(1, length(scleria_probocc_2.8km), by=2)] %>% terra::rast()
# set & simplify names
names(scleria_probocc_2.8km) <- varnames(scleria_probocc_2.8km)
names(scleria_probocc_2.8km) <- sub(' ', '_', names(scleria_probocc_2.8km) %>% stringr::word(1, 2))
varnames(scleria_probocc_2.8km) <- sub(' ', '_', varnames(scleria_probocc_2.8km) %>% stringr::word(1, 2))
# save stack
writeRaster(scleria_probocc_2.8km, 'results/maps/scleria_probocc_2.8km.tiff', overwrite=TRUE)

# get richness
scleria_richness_2.8km <- terra::app(scleria_probocc_2.8km, fun=sum, na.rm=TRUE)
# change raster names
names(scleria_richness_2.8km) <- 'scleria_richness'
# maintain NAs
scleria_richness_2.8km <- sum(scleria_richness_2.8km, r, na.rm=FALSE)
# save
writeRaster(scleria_richness_2.8km, 'results/maps/scleria_richness_2.8km.tiff', overwrite=TRUE)



# AGGREGATED

# import extended SDMs
scleria_probocc_50km <- list.files('SDM/results/eqearth_predictions/50km_agg/', full.names=TRUE)
# exclude .aux files and prepare stack
scleria_probocc_50km <- scleria_probocc_50km[seq(1, length(scleria_probocc_50km), by=2)] %>% terra::rast()
# set & simplify names
names(scleria_probocc_50km) <- varnames(scleria_probocc_50km)
names(scleria_probocc_50km) <- sub(' ', '_', names(scleria_probocc_50km) %>% stringr::word(1, 2))
varnames(scleria_probocc_50km) <- sub(' ', '_', varnames(scleria_probocc_50km) %>% stringr::word(1, 2))
# save stack
writeRaster(scleria_probocc_50km, 'results/maps/scleria_probocc_50km.tiff', overwrite=TRUE)

# get richness
scleria_richness_50km <- terra::app(scleria_probocc_50km, fun=sum, na.rm=TRUE)
# change raster names
names(scleria_richness_50km) <- 'scleria_richness'
# maintain NAs
scleria_richness_50km <- sum(scleria_richness_50km, r_agg, na.rm=FALSE)
# save
writeRaster(scleria_richness_50km, 'results/maps/scleria_richness_50km.tiff', overwrite=TRUE)



# COMPARE RESULTS TO PREVIOUS SDMs (MAXENT with less variables) ####

# previous analysis
temp <- terra::rast('results/maps/old_2023/exp_alpha_div_rasters.tiff') 
temp <- temp$richness
crs(temp) <- '+proj=eqearth'

# resample to resolution of new results and data.frame
temp <- temp %>% terra::resample(scleria_richness_50km) %>% as.data.frame(xy=TRUE)

# extract values of new richness raster
temp$richness_new <- terra::extract(scleria_richness_50km, temp[,c('x','y')]) %>% deframe()

# plot
ggplot(aes(y=richness_new, x=richness), data=temp) +
  geom_point() +
  xlab('Richness MAXENT 2023') + ylab('Richness STACKED 2025') +
  geom_smooth(method='lm') +
  theme_bw()

# there is a very good positive correlation and range is the same but, as we are not using thresholded values, results are more conservative now


