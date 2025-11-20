

# Lets explore the distribution of Scleria subgenera

source('scripts/1) import and prepare data.R')
library(ggpubr)


# import data
scleria_probocc_50km <- terra::rast('results/maps/scleria_probocc_50km.tiff')
names(scleria_probocc_50km) <- sub(".*_", "", names(scleria_probocc_50km))
names(scleria_probocc_50km)[names(scleria_probocc_50km)=='baroni-clarkei'] <- 'baroniclarkei'
names(scleria_probocc_50km)[names(scleria_probocc_50km)=='flagellum-nigrorum'] <- 'flagellumnigrorum'
names(scleria_probocc_50km)[names(scleria_probocc_50km)=='novae-hollandiae'] <- 'novaehollandiae'



# map subgenera
mapbrow <- scleria_probocc_50km %>% terra::subset(scl_taxonomy$epithet[scl_taxonomy$subgenus=='Browniae']) %>%
  terra::app(fun=sum, na.rm=TRUE)
maptrac <- scleria_probocc_50km %>% terra::subset(scl_taxonomy$epithet[scl_taxonomy$subgenus=='Trachylomia']) %>%
  terra::app(fun=sum, na.rm=TRUE)
mapscle <- scleria_probocc_50km %>% terra::subset(scl_taxonomy$epithet[scl_taxonomy$subgenus=='Scleria']) %>%
  terra::app(fun=sum, na.rm=TRUE)
maphypo <- scleria_probocc_50km %>% terra::subset(scl_taxonomy$epithet[scl_taxonomy$subgenus=='Hypoporum']) %>%
  terra::app(fun=sum, na.rm=TRUE)

# map
par(mfrow=c(2,2))
plot(mapbrow, main='subg. Browniae'); lines(world_lines)
plot(maptrac, main='subg. Trachylomia'); lines(world_lines)
plot(mapscle, main='subg. Scleria'); lines(world_lines)
plot(maphypo, main='subg. Hypoporum'); lines(world_lines)



# map sections
temp <- unique(scl_taxonomy$section)

# map
par(mfrow=c(3,2))
for (s in 1:length(temp)) {
  plot(scleria_probocc_50km %>% terra::subset(scl_taxonomy$epithet[scl_taxonomy$section==temp[s]]) %>%
  terra::app(fun=sum, na.rm=TRUE), main=paste('section', temp[s], sep=' '))
}

