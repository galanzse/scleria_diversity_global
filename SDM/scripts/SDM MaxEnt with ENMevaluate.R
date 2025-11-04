

# Run MaxEnt for species with 15 to 49 observations
source('SDM/scripts/SDM organize analyses.R')

library(ENMeval)
library(ks) # kde
library(spatialEco) # raster.transformation
# library(ecospat)
# library(remotes)
# install_github("r-spatial/sf")
# library(compiler, lib.loc = "C:/Program Files/R/R-4.5.0/library")
# library(tools, lib.loc = "C:/Program Files/R/R-4.5.0/library")



getwd()
# setwd('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity')


# define species and occurrences
scl_occ_max <- scl_sdm_alg %>% subset(algorithm=='maxent')


# loop to model distributions and obtain sites with high habitat suitability for posterior analyses
allmodels <- list()
bestmodels <- list()
varimp_bestmodels <- list()

# plot parameters
cuts <- seq(0, 1, 0.1)
pal <- colorRampPalette(c("white","black"))


for (s1 in 1:nrow(scl_occ_max)) {

  # occurrences
  occ_temp <- scl_occurrences %>% filter(scientific_name==scl_occ_max$species[s1])

  # points
  pts_p <- vect(occ_temp, geom=c('x','y'), 'epsg:4326')

  # model in ecoregions that intersect its EOO
  aoi_eco <- AOI %>% terra::mask(terra::convHull(pts_p)) # select geometries of x that intersect with the geometries of y

  # crop predictors
  r <- crop_unc_rst_predictors_2.5min %>% terra::crop(aoi_eco, mask=T)
  

  # plot(r$drainage_density)
  # lines(aoi_eco, col='red')
  # lines(countries)
  # points(pts_p, col='blue')
  
  
  # Bias file: Target-group background (https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13442)
  target_points <- data_final[['occurrences']][,c('x','y')] %>% unique() %>% 
    terra::vect(geom=c('x','y'), 'epsg:4326') %>%
    terra::intersect(aoi_eco) # all Scleria observations for our study area
  
  # species modeled with ensemble have large distribution ranges so we need to restrict the generation of bg points to sites not too far to presences (e.g. 500km) so bg points are sites environmentally kind of close to the presence points
  bff_p <- pts_p %>% terra::buffer(width=500000) %>% terra::aggregate(dissolve=TRUE) # 500km
  target_points <- target_points %>% terra::crop(bff_p)
  
  # 2d kernel density estimation 
  target_density <- ks::kde(geom(target_points)[,c('x','y')]) 
  
  # resample + crop to ground
  target_raster <- raster::raster(target_density) %>% terra::rast() %>%
    terra::resample(r$drainage_density, method='bilinear') %>%
    terra::crop(bff_p, mask=T)

  # normalize
  target_raster <- target_raster - min(target_raster[], na.rm=T)
  target_raster <- raster.transformation(target_raster, trans="norm") 

  # select a maximum of 1000 background points
  bg <- as.data.frame(target_raster, xy=T) %>% subset(layer>0.1)
  if (nrow(bg) > 1000) { bg <- bg[sample(1:nrow(bg), prob=bg$layer, size=1000, replace=F),] }

  # plot(target_raster)
  # points(bg[,c('x','y')], col='orange')
  # points(target_points, col='red', pch=16)
  # points(pts_p, col='blue', pch=16)

  
  # Get best model: ENMeval
  modeval <- ENMeval::ENMevaluate(occs = occ_temp[,c('x','y')], 
                         envs = r,
                         bg= bg[,c('x','y')],
                         algorithm = 'maxent.jar',
                         tune.args = list(fc = c('L', 'Q', 'LQ', 'LQH'), rm = 1:5), 
                         partitions = "randomkfold", partition.settings = list(kfolds = 5), 
                         doClamp = TRUE)

  
  # get results and add number of features
  modeval_results <- modeval@results
  modeval_results$nfeat <- nchar(as.character(modeval_results$fc)) # number of features

  # get best model (lowest AICc > lowest ncoef > lowes nfeat > any)
  bestmodel <- modeval_results %>% subset(delta.AICc<2)
  if (nrow(bestmodel)>1) { bestmodel <- bestmodel %>% subset(ncoef==min(bestmodel$ncoef, na.rm=T)) }
  if (nrow(bestmodel)>1) { bestmodel <- bestmodel %>% subset(nfeat==min(bestmodel$nfeat, na.rm=T)) }
  if (nrow(bestmodel)>1) { bestmodel <- bestmodel[1,] }
  
  # best model index
  best_mod_ind <- as.numeric(rownames(bestmodel))
  
  # variable importance
  varimp_bm <- eval.variable.importance(modeval)[[best_mod_ind]]

  # prediction best model
  proj_temp <- modeval@predictions[[best_mod_ind]]

  # plot(proj_temp, main="Model with lowest AICc")
  # points(pts_p, col='blue')
  

  # SAVE OUTPUT
  
  # prediction
  writeRaster(proj_temp, paste('SDM/results/predictions/', scl_occ_max$species[s1], '.tiff', sep=''), overwrite=T)
  
  # map
  png(paste('SDM/results/maps/', scl_occ_max$species[s1], '.png', sep=''), bg="white")
  plot(proj_temp, breaks=cuts, col=pal(10), main=scl_occ_max$species[s1])
  lines(aoi_eco, col='black')
  points(pts_p, col='green', pch=16)
  dev.off()

  # model output:
    # all models
  modeval_results$species <- scl_occ_max$species[s1] # add species name
  allmodels[[s1]] <- modeval_results # save as indexed in loop
  save(allmodels, file='SDM/results/maxent/allmodels.Rdata') # overwrite results list
    # best model
  bestmodel$species <- scl_occ_max$species[s1]
  bestmodels[[s1]] <- bestmodel
  save(bestmodels, file='SDM/results/maxent/bestmodels.Rdata')
    # variable importance
  varimp_bm$species <- scl_occ_max$species[s1]
  varimp_bestmodels[[s1]] <- varimp_bm
  save(varimp_bestmodels, file='SDM/results/maxent/varimp_bestmodels.Rdata')

  
  # progress
  print(paste(as.character(scl_occ_max$species[s1]), '--- Done!   ->   ', round(s1/nrow(scl_occ_max)*100, 2), '%'))

}


