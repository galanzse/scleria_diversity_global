

# Run IWD for species with 4 to 9 observations
source('SDM/scripts/SDM organize analyses.R')


# 1/ Selection of background points - sampling bias (target-group approach)
# 2/ Fit model and retain best parametrization
# 3/ Save presence probability

library(dismo) # geoIDW
library(ks) # kde
library(spatialEco) # raster.transformation
library(gstat) # required for geoIDW
library(raster) # raster
library(pROC) # roc


# results dataframe
scl_occ_idw <- scl_sdm_alg %>% subset(algorithm=='IDW')
scl_occ_idw$AUC <- NA

# plot parameters
cuts <- seq(0, 1, 0.1)
pal <- colorRampPalette(c("white","black"))

# loop
for (i in 1:nrow(scl_occ_idw)) {
  
  # presence df and vect
  df_p <- scl_occurrences %>% subset(scientific_name==scl_occ_idw$species[i]) %>% dplyr::select(x,y)
  pts_p <- df_p %>% terra::vect(geom=c('x','y'), 'epsg:4326')

  # model in ecoregions that intersect the species EOO
  aoi_eco <- AOI %>% terra::mask(terra::convHull(pts_p)) %>% terra::aggregate(dissolve=TRUE)
  
  # base raster
  r <- crop_unc_rst_predictors_2.5min['elevation'] %>%
    terra::crop(aoi_eco, mask=T) %>% raster::raster()
  
  # bias file: target-group background (https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13442)
  target_points <- scl_occurrences %>% unique() %>%
    terra::vect(geom=c('x','y'), 'epsg:4326') %>%
    terra::intersect(aoi_eco)
  
  # 2d kernel density estimation
  target_density <- ks::kde(geom(target_points)[,c('x','y')])
  
  # resample + crop
  target_raster <- raster::raster(target_density) %>% terra::rast() %>%
    terra::resample(terra::rast(r), method='bilinear') %>%
    terra::crop(aoi_eco, mask=T)
  
  # normalize
  target_raster <- target_raster - min(target_raster[], na.rm=T)
  target_raster <- spatialEco::raster.transformation(target_raster, trans="norm")
  
  # generate background points to achieve a minimum prevalence of 10% (i.e., 1:10 ratio between the number of cells occupied by species records and the number of non-occupied cells in the study ecoregion)
  bg <- as.data.frame(target_raster, xy=T) %>% subset(layer>0.1)
  colnames(bg)[colnames(bg)=='layer'] <- 'prob'
  
  # run model and evaluate using leave-one-out
  n <- nrow(df_p)
  predicted_values <- numeric(n)
  
  # save each raster prediction in stack to average predictions
  # use r as template
  prd <- terra::rast(r)
  names(prd) <- 'template'
  
  # create as many predictions as allowed by the data
  for (i2 in 1:n) {
    # Leave-one-out: remove i-th point
    train_data <- df_p[-i2, ]
    test_point <- df_p[i2, ]
    # select background points
    bg1 <- bg[sample(x=1:nrow(bg), prob=bg$prob, size=nrow(df_p)*10, replace=F),][,c('x','y')]
    # Fit geoIDW model
    model <- dismo::geoIDW(train_data, bg1)
    # Predict at the left-out point
    pred1 <- predict(model, test_point)
    predicted_values[i2] <- pred1
    # predict using last model and combine predictions
    prd <- prd %>% c(predict(model, r) %>% terra::rast() %>% terra::crop(aoi_eco, mask=T))
  }

  # remove first 'template' layer
  prd$template <- NULL
  
  # average predictions for posterior analyses
  prd <- mean(prd)
  
  # save
  writeRaster(prd, paste('SDM/results/predictions/', scl_occ_idw$species[i], '.tiff', sep=''), overwrite=T)
  
  # combine observed (1 for presence) and predicted values
  observed <- rep(1, n) # presence-only
  results <- data.frame(obs=observed, pred=predicted_values)
  
  # evaluate using a different set of background points (if not AUC always equals 1)
  bg2 <- bg[sample(x=1:nrow(bg), prob=bg$prob, size=nrow(df_p)*10, replace=F),][,c('x','y')]
  bg_pred <- extract(prd, bg2, ID=F) %>% deframe()
  scl_occ_idw$AUC[i] <- pROC::roc(c(rep(1, n), rep(0, nrow(bg2))), c(results$pred, bg_pred))$auc
  
  # map
  png(paste('SDM/results/maps/', scl_occ_idw$species[i], '.png', sep=''), bg="white")
  plot(prd, breaks=cuts, col=pal(10), main=scl_occ_idw$species[i])
  lines(aoi_eco, col='black')
  points(df_p, col='green', pch=16)
  dev.off()
  
  print(paste(scl_occ_idw$species[i], '--- Done!', '[', round(i/nrow(scl_occ_idw),2)*100, '% ]'))
  
  # save results
  write.table(scl_occ_idw, 'SDM/results/IDW_AUC.txt')
  
}


