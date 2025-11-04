

# Identify uncorrelated variables


library(tidyverse)
library(terra)
library(corrplot) # corrplot
library(usdm) # vif


# data
rst_predictors_2.5min <- terra::rast('SDM/rst_predictors_2.5min.tiff')
names(rst_predictors_2.5min)

# remove categorical raster because we are modelling within ecoregions (already considered)
# and landcover is included in continuous rasters such as treecover and herbaceous_vegetation
rst_predictors_2.5min$biomes <- NULL
rst_predictors_2.5min$ecoregions <- NULL
rst_predictors_2.5min$landcover <- NULL


# extract values of stack
load("C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/data/data_final.RData")
pts_scleria <- data_final$occurrences %>% terra::vect(geom=c('x','y'), crs='epsg:4326')
temp <- rst_predictors_2.5min %>% terra::extract(pts_scleria, cells=F, xy=F, ID=F)

# remove NAs
temp <- na.omit(temp)


# calculate correlation matrix
var.cor <- temp %>% cor(method='spearman', use='pairwise.complete.obs')
  
# graphical display of a correlation matrix
corrplot(var.cor, type="upper", method="number", tl.cex=1, cl.cex=1, cl.ratio=0.1, col=COL2('RdBu', 10))

# matrix to dataframe
cor.df <- as.data.frame(var.cor)

# keep upper tri
lower <- var.cor
lower[lower.tri(var.cor, diag=TRUE)] <- ""
lower.df <- as.data.frame(lower)

# correlation matrix to distance matrix
var.dist <- abs(as.dist(cor.df))

# build dendrogram: distance is inversely proportional to the value of the correlation coefficient (shorter distance = higher correlation)
var.cluster <- hclust(1-var.dist)

# plot
plot(var.cluster)
abline(h=1-0.8, lty=2, lwd=2, col="red")


# remove highly correlated variables
rst_predictors_2.5min$bio_3 <- NULL
rst_predictors_2.5min$bio_4 <- NULL
rst_predictors_2.5min$bio_13 <- NULL
rst_predictors_2.5min$bio_6 <- NULL
rst_predictors_2.5min$bio_9 <- NULL
rst_predictors_2.5min$bio_11 <- NULL
rst_predictors_2.5min$bio_14 <- NULL
rst_predictors_2.5min$bio_15 <- NULL
rst_predictors_2.5min$pop_density <- NULL
rst_predictors_2.5min$bio_5 <- NULL
rst_predictors_2.5min$slope <- NULL
rst_predictors_2.5min$aridity <- NULL
rst_predictors_2.5min$water_content_10 <- NULL


# remove other variables with little impact or expected causality
names(rst_predictors_2.5min)
rst_predictors_2.5min$aspect <- NULL
rst_predictors_2.5min$evapotrans <- NULL
rst_predictors_2.5min$water_content_1500 <- NULL # 33 is the field capacity
rst_predictors_2.5min$bio_10 <- NULL
rst_predictors_2.5min$bio_16 <- NULL
rst_predictors_2.5min$bio_19 <- NULL
rst_predictors_2.5min$flooded_vegetation <- NULL


# redo dataframe of extracted values
temp <- rst_predictors_2.5min %>% terra::extract(pts_scleria, cells=F, xy=F, ID=F) %>% na.omit()

# run vifstep on final 27 variables for final selection
# https://quantifyinghealth.com/vif-threshold/
vif_temp <- usdm::vifstep(temp, th=2.5)
vif_temp@results$Variables

# save uncorrelated predictors
# unc_rst_predictors_2.5min <- rst_predictors_2.5min %>% subset(vif_temp@results$Variables)
# names(unc_rst_predictors_2.5min)
# writeRaster(unc_rst_predictors_2.5min, 'SDM/unc_rst_predictors_2.5min.tiff', overwrite=TRUE)


