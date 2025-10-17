

# CALCULATE RICHNESS and MORAN'S I FOR THE 50X50KM DATASET
# PREVALENCE IS 100% FOR THE NEW DATASET

source('scripts/1) import and prepare data.R')
source('functions/fun_autocor.R')

library(ggpubr)
library(ape) # moran.i
library(geosphere)


# MAT as base raster
MAT_ee <- terra::rast('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_1.tif') %>% terra::project('+proj=eqearth')
cellSize(MAT_ee, unit='km') %>% sqrt() # cell side length

# aggregate to 50*50km (2500km2)
MAT_ee  %>% terra::aggregate(18) %>% cellSize(unit='km') %>% sqrt()
MAT_ee <- MAT_ee %>% terra::aggregate(18)

# back to lanlon
MAT_ll <- terra::rast('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/worldclim/wc2.1_2.5m historical_19702000/wc2.1_2.5m_bio_1.tif') %>% terra::aggregate(18)


# occurrences - dataframe of species occurrences (y/x) 
pts_scleria_ll <- expected_assemblages %>%
  pivot_longer(4:ncol(expected_assemblages), names_to='species', values_to='abundance') %>%
  subset(abundance>0) %>%
  terra::vect(geom=c('x','y'), '+proj=eqearth') %>%
  terra::project('epsg:4326')

# project occurrences same as raster
pts_scleria_ee <- pts_scleria_ll %>% terra::project('+proj=eqearth')


# extract cell ID occupied per species
df_pts_scleria_ee <- MAT_ee %>% terra::extract(pts_scleria_ee, cells=T, xy=T) %>% dplyr::select(cell, x, y)
df_pts_scleria_ee$species <- pts_scleria_ee$species
  
# calculate mean species richness
df_pts_scleria_ee %>% group_by(cell) %>% summarise(temp=n_distinct(species)) %>% dplyr::select(temp) %>% deframe() %>% mean()
df_pts_scleria_ee %>% group_by(cell) %>% summarise(temp=n_distinct(species)) %>% dplyr::select(temp) %>% deframe() %>% sd()



# calculate Moran's I
# https://cran.r-project.org/web/packages/ape/vignettes/MoranI.pdf

df_pts_scleria_ll <- pts_scleria_ll
df_pts_scleria_ee <- MAT_ee %>% terra::extract(pts_scleria_ee, cells=T, xy=T) %>% dplyr::select(cell, x, y)
df_pts_scleria_ee$species <- pts_scleria_ee$species


# retain unique cells
pts_scleria_ll$cell <- NULL; pts_scleria_ll$species <- NULL; pts_scleria_ll$abundance <- NULL
pts_scleria_ll <- unique(pts_scleria_ll)

# subset some points to it does not crack
temp <- pts_scleria_ll[sample(1:nrow(pts_scleria_ll),500),]

# extract values from the predictor variable for the points of presence
xy.variables <- na.omit(
  data.frame(geom(temp)[,c('x','y')],
             terra::extract(x=MAT_ll, y=temp, df=TRUE, cells=FALSE, ID=FALSE)
  )
)
  
# rename variable to call it later
colnames(xy.variables)[!(colnames(xy.variables)%in%c('x','y'))] <- 'var'
  
# matrix of inverse distances between occurrences
xy.distancias <- 1/ (geosphere::distm(x=xy.variables[, c("x","y")], fun=distGeo) / 1000)
xy.distancias[!is.finite(xy.distancias)] <- 0 # replace Inf with zero
diag(xy.distancias) <- 0 # diagonal to 0
  
# Moran's I 
ape::Moran.I(xy.variables$var, xy.distancias, na.rm=TRUE)


