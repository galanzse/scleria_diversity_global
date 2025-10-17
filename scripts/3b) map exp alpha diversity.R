

# MAP ALPHADIVERSITY FOR EXPECTED ASSEMBLAGES


source('scripts/1) import and prepare data.R')

library(SpatialPack) # modified.ttest


# data
exp_alpha_div <- read.csv("results/exp_alpha_div.txt", sep="")
world_grid



# dataframe to map #### 

# convert n_ann and n_per into an index
exp_alpha_div$lifeform_index <- c(exp_alpha_div$n_per - exp_alpha_div$n_ann) / c(exp_alpha_div$n_per + exp_alpha_div$n_ann)



# template to fill
alpha_div_rasters <- rep(world_grid, 14)
names(alpha_div_rasters) <- c("richness", "shannon", "cwm_height", "cwm_blade", "cwm_nutlet", "lifeform_index", "Frich", "Fmpd", "Prich","Pmpd","perc_Frich", "perc_Fmpd", "perc_Prich", "perc_Pmpd")

# fill values
alpha_div_rasters$richness <- terra::rast(exp_alpha_div[,c('x','y','richness')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$shannon <- terra::rast(exp_alpha_div[,c('x','y','shannon')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)

alpha_div_rasters$cwm_height <- terra::rast(exp_alpha_div[,c('x','y','cwm_height')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$cwm_blade <- terra::rast(exp_alpha_div[,c('x','y','cwm_blade')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$cwm_nutlet <- terra::rast(exp_alpha_div[,c('x','y','cwm_nutlet')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$lifeform_index <- terra::rast(exp_alpha_div[,c('x','y','lifeform_index')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)

alpha_div_rasters$Frich <- terra::rast(exp_alpha_div[,c('x','y','Frich')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$Fmpd <- terra::rast(exp_alpha_div[,c('x','y','Fmpd')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$Prich <- terra::rast(exp_alpha_div[,c('x','y','Prich')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$Pmpd <- terra::rast(exp_alpha_div[,c('x','y','Pmpd')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)

alpha_div_rasters$perc_Frich <- terra::rast(exp_alpha_div[,c('x','y','perc_Frich')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$perc_Fmpd <- terra::rast(exp_alpha_div[,c('x','y','perc_Fmpd')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$perc_Prich <- terra::rast(exp_alpha_div[,c('x','y','perc_Prich')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)
alpha_div_rasters$perc_Pmpd <- terra::rast(exp_alpha_div[,c('x','y','perc_Pmpd')], crs='+proj=eqearth') %>%
  terra::extend(world_grid)

# save rasters
writeRaster(alpha_div_rasters, 'results/maps/alpha_div_expected_50/exp_alpha_div_rasters.tiff', overwrite=TRUE)



# prepare maps for representation ####

for (s in 1:dim(alpha_div_rasters)[3]) {
  
  pdf(file = paste("results/maps/alpha_div_expected_50/",names(alpha_div_rasters)[s],".pdf",sep=""),
      width = 9.30, # The width of the plot in inches
      height = 5.74) # The height of the plot in inches
  alpha_div_rasters[[s]] %>%
    plot(main=names(alpha_div_rasters)[s], col=rev(grDevices::heat.colors(20)[1:15]))
  lines(world_lines)
  
  dev.off()
  
}


# SES
par(mfrow=c(2,2))

temp <- alpha_div_rasters$perc_Frich
temp[temp>=0.95] <- 1; temp[temp<=0.05] <- 2; temp[temp<0.95 & temp>0.05] <- 3
plot(temp, col=c('indianred','blue','snow3'), legend=NULL, main='Significant Frich')
lines(world_lines, alpha=0.5)

temp <- alpha_div_rasters$perc_Fmpd
temp[temp>=0.95] <- 1; temp[temp<=0.05] <- 2; temp[temp<0.95 & temp>0.05] <- 3
plot(temp, col=c('indianred','blue','snow3'), legend=NULL, main='Significant Fmpd')
lines(world_lines, alpha=0.5)

temp <- alpha_div_rasters$perc_Prich
temp[temp>=0.95] <- 1; temp[temp<=0.05] <- 2; temp[temp<0.95 & temp>0.05] <- 3
plot(temp, col=c('indianred','blue','snow3'), legend=NULL, main='Significant Prich')
lines(world_lines, alpha=0.5)

temp <- alpha_div_rasters$perc_Pmpd
temp[temp>=0.95] <- 1; temp[temp<=0.05] <- 2; temp[temp<0.95 & temp>0.05] <- 3
plot(temp, col=c('indianred','blue','snow3'), legend=NULL, main='Significant Pmpd')
lines(world_lines, alpha=0.5)


