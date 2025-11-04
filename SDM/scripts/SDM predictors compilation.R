

# Prepare raster stack for SMs: resample at 2.5m (20km2 at equator)


library(tidyverse)
library(terra)
library(dismo)
# library(rnaturalearthdata)
# library(ecospat)


# import occurrence data (already curated)
load("C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/data/data_final.RData")


# define region of interest (roi)
scl_occurrences <- data_final$occurrences
scl_points <- terra::vect(scl_occurrences, geom=c('x','y'), crs='epsg:4326')
scl_hull <- scl_points %>% terra::hull(type="concave_ratio", param=0.5) %>% terra::buffer(width=500000) # 500km

plot(scl_hull, col='red'); points(scl_points)
# 
# 
# 
# import rasters, resample using wc2.1_2.5min_elev and export to world_rasters_rs2.5
# elevation <- terra::rast('C:/Users/javie/Desktop/world_rasters/elevation/wc2.1_2.5m_elev.tif') %>% terra::crop(scl_hull, mask=T)
# names(elevation) <- 'elevation'
# writeRaster(elevation, 'C:/Users/javie/Desktop/world_rasters_rs2.5/elevation.tiff', overwrite=TRUE)
# 
# 
# aspect <- terra::rast('C:/Users/javie/Desktop/world_rasters/elevation/aspect_2.5m.tif') %>% terra::crop(scl_hull, mask=T)
# names(aspect) <- 'aspect'
# writeRaster(aspect, 'C:/Users/javie/Desktop/world_rasters_rs2.5/aspect.tiff', overwrite=TRUE)
# 
# 
# ruggedness <- terra::rast('C:/Users/javie/Desktop/world_rasters/elevation/ruggedness_2.5m.tif') %>% terra::crop(scl_hull, mask=T)
# names(ruggedness) <- 'ruggedness'
# writeRaster(ruggedness, 'C:/Users/javie/Desktop/world_rasters_rs2.5/ruggedness.tiff', overwrite=TRUE)
# 
# 
# slope <- terra::rast('C:/Users/javie/Desktop/world_rasters/elevation/slope_2.5m.tif') %>% terra::crop(scl_hull, mask=T)
# names(slope) <- 'slope'
# writeRaster(slope, 'C:/Users/javie/Desktop/world_rasters_rs2.5/slope.tiff', overwrite=TRUE)
# 
# 
# aridity <- terra::rast('C:/Users/javie/Desktop/world_rasters/aridity/ai_v3_yr.tif') %>% terra::crop(scl_hull, mask=T)
# names(aridity) <- 'aridity'
# aridity <- aridity %>% terra::resample(elevation, method='bilinear')
# writeRaster(aridity, 'C:/Users/javie/Desktop/world_rasters_rs2.5/aridity.tiff', overwrite=TRUE)
# 
# 
# evapotrans <- terra::rast('C:/Users/javie/Desktop/world_rasters/aridity/et0_v3_yr.tif') %>% terra::crop(scl_hull, mask=T)
# names(evapotrans) <- 'evapotrans'
# evapotrans <- evapotrans %>% terra::resample(elevation, method='bilinear')
# writeRaster(evapotrans, 'C:/Users/javie/Desktop/world_rasters_rs2.5/evapotrans.tiff', overwrite=TRUE)
# 
# 
# ecoregions <- terra::rast('C:/Users/javie/Desktop/world_rasters/ecoregions/ecoregions2017/rst_ecoregions.tif') %>% terra::crop(scl_hull, mask=T)
# names(ecoregions) <- 'ecoregions'
# ecoregions <- ecoregions %>% terra::resample(elevation, method='mode')
# writeRaster(ecoregions, 'C:/Users/javie/Desktop/world_rasters_rs2.5/ecoregions.tiff', overwrite=TRUE)
# 
# 
# pop_density <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/population_density/Global_2000_PopulationDensity30sec_GPWv4.tiff",
#          "C:/Users/javie/Desktop/world_rasters/population_density/Global_2005_PopulationDensity30sec_GPWv4.tiff",
#          "C:/Users/javie/Desktop/world_rasters/population_density/Global_2010_PopulationDensity30sec_GPWv4.tiff",
#          "C:/Users/javie/Desktop/world_rasters/population_density/Global_2015_PopulationDensity30sec_GPWv4.tiff",
#          "C:/Users/javie/Desktop/world_rasters/population_density/Global_2020_PopulationDensity30sec_GPWv4.tiff")) %>%
#   terra::crop(scl_hull, mask=T)
# pop_density <- pop_density %>% terra::app(mean, na.rm=TRUE) %>% terra::resample(elevation, method='bilinear')
# names(pop_density) <- 'pop_density'
# writeRaster(pop_density, 'C:/Users/javie/Desktop/world_rasters_rs2.5/pop_density.tiff', overwrite=TRUE)
# 
# 
# treecover <- terra::rast('C:/Users/javie/Desktop/world_rasters/percent_tree_cover/gm_ve_v2.tif') %>%
#   terra::crop(scl_hull, mask=T)
# names(treecover) <- 'treecover'
# treecover[treecover>100] <- NA # water bodies and no data
# treecover <- treecover %>% resample(elevation, method='bilinear')
# writeRaster(treecover, 'C:/Users/javie/Desktop/world_rasters_rs2.5/treecover.tiff', overwrite=TRUE)
# 
# 
# crs1 <- terra::rast("C:/Users/javie/Desktop/world_rasters/human_footprint/hfp2000_merisINT.tif") %>% crs()
# human_footprint <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/human_footprint/hfp2000_merisINT.tif",
#   "C:/Users/javie/Desktop/world_rasters/human_footprint/hfp2005_merisINT.tif",
#   "C:/Users/javie/Desktop/world_rasters/human_footprint/hfp2010_merisINT.tif",
#   "C:/Users/javie/Desktop/world_rasters/human_footprint/hfp2013_merisINT.tif")) %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T)
# human_footprint <- human_footprint %>% terra::app(mean, na.rm=TRUE) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(human_footprint) <- 'human_footprint'
# writeRaster(human_footprint, 'C:/Users/javie/Desktop/world_rasters_rs2.5/human_footprint.tiff', overwrite=TRUE)
# 
# 
# biomes <- vect('C:/Users/javie/Desktop/world_rasters/ecoregions/ecoregions2017/Ecoregions2017.shp')['BIOME_NUM']
# names(biomes) <- 'biomes'
# biomes <- terra::rasterize(x=biomes, y=elevation, field='biomes', fun='min')
# biomes <- biomes %>% terra::crop(scl_hull, mask=T)
# writeRaster(biomes, 'C:/Users/javie/Desktop/world_rasters_rs2.5/biomes.tiff', overwrite=TRUE)
# 
# 
# drainage_density <- terra::rast('C:/Users/javie/Desktop/world_rasters/global drainage density/LCS_Dd_global.tif') %>% terra::crop(scl_hull, mask=T)
# names(drainage_density) <- 'drainage_density'
# drainage_density <- drainage_density %>% terra::resample(elevation, method='bilinear')
# writeRaster(drainage_density, 'C:/Users/javie/Desktop/world_rasters_rs2.5/drainage_density.tiff', overwrite=TRUE)
# 
# 
# crs1 <- terra::rast("C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_0-5cm_mean_5000.tif") %>% crs()
# soil_nitrogen <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_0-5cm_mean_5000.tif",
#                         "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_5-15cm_mean_5000.tif",
#                         "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_15-30cm_mean_5000.tif",
#                         "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_30-60cm_mean_5000.tif",
#                         "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_60-100cm_mean_5000.tif",
#                         "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/nitrogen/nitrogen_100-200cm_mean_5000.tif")) %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T)
# soil_nitrogen <- soil_nitrogen %>% terra::app(sum, na.rm=TRUE) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(soil_nitrogen) <- 'soil_nitrogen'
# writeRaster(soil_nitrogen, 'C:/Users/javie/Desktop/world_rasters_rs2.5/soil_nitrogen.tiff', overwrite=TRUE)
# 
# 
# organic_c_stock <- terra::rast('C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/ocs/ocs_0-30cm_mean_5000.tif') %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(organic_c_stock) <- 'organic_c_stock'
# writeRaster(organic_c_stock, 'C:/Users/javie/Desktop/world_rasters_rs2.5/organic_c_stock.tiff', overwrite=TRUE)
# 
# 
# soil_ph <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/phh2o/phh2o_0-5cm_mean_5000.tif",
#                   "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/phh2o/phh2o_5-15cm_mean_5000.tif",
#                   "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/phh2o/phh2o_15-30cm_mean_5000.tif",
#                   "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/phh2o/phh2o_30-60cm_mean_5000.tif",
#                   "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/phh2o/phh2o_60-100cm_mean_5000.tif",
#                   "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/phh2o/phh2o_100-200cm_mean_5000.tif")) %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T)
# soil_ph <- soil_ph %>% terra::app(mean, na.rm=TRUE) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(soil_ph) <- 'soil_ph'
# writeRaster(soil_ph, 'C:/Users/javie/Desktop/world_rasters_rs2.5/soil_ph.tiff', overwrite=TRUE)
# 
# 
# water_content_10 <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0010/wv0010_0-5cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0010/wv0010_5-15cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0010/wv0010_15-30cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0010/wv0010_30-60cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0010/wv0010_60-100cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0010/wv0010_100-200cm_mean_5000.tif")) %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T)
# water_content_10 <- water_content_10 %>% terra::app(sum, na.rm=TRUE) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(water_content_10) <- 'water_content_10'
# writeRaster(water_content_10, 'C:/Users/javie/Desktop/world_rasters_rs2.5/water_content_10.tiff', overwrite=TRUE)
# 
# 
# water_content_33 <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0033/wv0033_0-5cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0033/wv0033_5-15cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0033/wv0033_15-30cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0033/wv0033_30-60cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0033/wv0033_60-100cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv0033/wv0033_100-200cm_mean_5000.tif")) %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T)
# water_content_33 <- water_content_33 %>% terra::app(sum, na.rm=TRUE) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(water_content_33) <- 'water_content_33'
# writeRaster(water_content_33, 'C:/Users/javie/Desktop/world_rasters_rs2.5/water_content_33.tiff', overwrite=TRUE)
# 
# 
# water_content_1500 <- terra::rast(c("C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv1500/wv1500_0-5cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv1500/wv1500_5-15cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv1500/wv1500_15-30cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv1500/wv1500_30-60cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv1500/wv1500_60-100cm_mean_5000.tif",
#                            "C:/Users/javie/Desktop/world_rasters/isric_soilgrids_5km/wv1500/wv1500_100-200cm_mean_5000.tif")) %>%
#   terra::crop(terra::project(scl_hull, crs1), mask=T)
# water_content_1500 <- water_content_1500 %>% terra::app(sum, na.rm=TRUE) %>%
#   terra::project('epsg:4326') %>%
#   terra::resample(elevation, method='bilinear')
# names(water_content_1500) <- 'water_content_1500'
# writeRaster(water_content_1500, 'C:/Users/javie/Desktop/world_rasters_rs2.5/water_content_1500.tiff', overwrite=TRUE)
# 
# 
# flooded_vegetation <- terra::rast('C:/Users/javie/Desktop/world_rasters/Global 1km Consensus Land Cover/consensus_full_class_8.tif') %>%
#   terra::crop(scl_hull, mask=T)
# names(flooded_vegetation) <- 'flooded_vegetation'
# flooded_vegetation <- flooded_vegetation %>% terra::resample(elevation, method='bilinear')
# writeRaster(flooded_vegetation, 'C:/Users/javie/Desktop/world_rasters_rs2.5/flooded_vegetation.tiff', overwrite=TRUE)
# 
# 
# herbaceous_vegetation <- terra::rast('C:/Users/javie/Desktop/world_rasters/Global 1km Consensus Land Cover/consensus_full_class_6.tif') %>%
#   terra::crop(scl_hull, mask=T)
# names(herbaceous_vegetation) <- 'herbaceous_vegetation'
# herbaceous_vegetation <- herbaceous_vegetation %>% terra::resample(elevation, method='bilinear')
# writeRaster(herbaceous_vegetation, 'C:/Users/javie/Desktop/world_rasters_rs2.5/herbaceous_vegetation.tiff', overwrite=TRUE)
# 
# 
# landcover <- terra::rast('C:/Users/javie/Desktop/world_rasters/percent_tree_cover/gm_lc_v3.tif') %>%
#   terra::crop(scl_hull, mask=T)
# names(landcover) <- 'landcover'
# landcover[landcover%in%c(19,20)] <- NA # ice, snow and water bodies
# landcover <- landcover %>% resample(elevation, method='min')
# writeRaster(landcover, 'C:/Users/javie/Desktop/world_rasters_rs2.5/landcover.tiff', overwrite=TRUE)
# 
# 
# collection effort in Scleria peaks in between 1970 and 2010, so lets use historical data
# worldclim1 <- rast(c(
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_19.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_1.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_2.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_3.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_4.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_5.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_6.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_7.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_8.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_9.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_10.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_11.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_12.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_13.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_14.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_15.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_16.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_17.tif",
#   "C:/Users/javie/Desktop/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_18.tif"
# )) %>%  terra::crop(scl_hull, mask=T)
# names(worldclim1) <- c("bio_19", "bio_1",  "bio_2",  "bio_3",  "bio_4", 
#                        "bio_5",  "bio_6",  "bio_7",  "bio_8",  "bio_9" ,
#                        "bio_10", "bio_11", "bio_12", "bio_13", "bio_14",
#                        "bio_15", "bio_16", "bio_17", "bio_18")
# writeRaster(worldclim1, 'C:/Users/javie/Desktop/world_rasters_rs2.5/worldclim.tiff')
# 
# 
# # merge into one file
# rst_predictors <- terra::rast(list.files('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/world_rasters_rs2.5', full.names=TRUE))
# names(rst_predictors)


# filter and save
# getwd()
# writeRaster(rst_predictors, 'C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/rst_predictors_2.5min.tiff', overwrite=T)


