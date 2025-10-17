

library(tidyverse)
library(terra)
library(SpatialPack) # modified.ttest


# obs alpha diversity ~ climate ####

# response variables 
alphadiv <- terra::rast('results/maps/alpha_div_observed_50/obs_alpha_div_rasters.tiff')[[1:9]] %>%
  terra::project('epsg:4326') %>%
  as.data.frame(xy=T)
names(alphadiv)

# bioclim
wc_bioclim <- terra::rast('SDM/rst_predictors_2.5min.tiff')[[c("bio_1","bio_7","bio_12","bio_15")]] %>%
  terra::aggregate(fact=12)
cellSize(wc_bioclim, unit='km') %>% sqrt() # ~50km2

# extract bioclim
pts_alphadiv <- terra::vect(alphadiv, geom=c('x','y'), 'epsg:4326')
alphadiv$bio_1 <- wc_bioclim$bio_1 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()
alphadiv$bio_7 <- wc_bioclim$bio_7 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()
alphadiv$bio_12 <- wc_bioclim$bio_12 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()
alphadiv$bio_15 <- wc_bioclim$bio_15 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()

test_alphadiv_bioclim <- data.frame(ALPHA=rep(c("richness","Frich","Fmpd","Prich","Pmpd",
                                                "cwm_height","cwm_blade","cwm_nutlet"), 5),
                                    BIOCLIM=c( rep('bio_1',8), rep('bio_7',8), rep('bio_12',8), rep('bio_15',8),
                                               rep('y',8)),
                                    F_=NA, ESS=NA, DF=NA, corr=NA, P_value=NA, B_=NA)

# transform latitude to reflect distance to the equator
alphadiv$y[alphadiv$y<0] <- alphadiv$y[alphadiv$y<0] * c(-1)

for (i in 1:nrow(test_alphadiv_bioclim)) {
  
  temp <- alphadiv[,c(test_alphadiv_bioclim$ALPHA[i], test_alphadiv_bioclim$BIOCLIM[i], 'x', 'y')] %>% na.omit()
  
  modtest1 <- modified.ttest(temp[,test_alphadiv_bioclim$ALPHA[i]],
                             temp[,test_alphadiv_bioclim$BIOCLIM[i]],
                             temp[,c('x','y')], nclass=13)
  
  test_alphadiv_bioclim$F_[i] <- modtest1$Fstat
  test_alphadiv_bioclim$ESS[i] <- modtest1$ESS
  test_alphadiv_bioclim$DF[i] <- modtest1$dof
  test_alphadiv_bioclim$corr[i] <- modtest1$corr
  test_alphadiv_bioclim$P_value[i] <- modtest1$p.value
  test_alphadiv_bioclim$B_[i] <- lm(temp[,1] ~ temp[,2])$coefficients[2]
  
  print(paste('--- ', round(i/nrow(test_alphadiv_bioclim),2)*100,'% ---', sep=''))
  
}

write.csv(test_alphadiv_bioclim, 'results/test_obs_alphadiv_bioclim.csv')



# exp alpha diversity ~ climate ####

# response variables 
alphadiv <- terra::rast('results/maps/alpha_div_expected_50/exp_alpha_div_rasters.tiff')[[1:9]] %>%
  terra::project('epsg:4326') %>%
  as.data.frame(xy=T)
names(alphadiv)

# bioclim
wc_bioclim <- terra::rast('SDM/rst_predictors_2.5min.tiff')[[c("bio_1","bio_7","bio_12","bio_15")]] %>%
  terra::aggregate(fact=12)
cellSize(wc_bioclim, unit='km') %>% sqrt() # ~50km2

# extract bioclim
pts_alphadiv <- terra::vect(alphadiv, geom=c('x','y'), 'epsg:4326')
alphadiv$bio_1 <- wc_bioclim$bio_1 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()
alphadiv$bio_7 <- wc_bioclim$bio_7 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()
alphadiv$bio_12 <- wc_bioclim$bio_12 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()
alphadiv$bio_15 <- wc_bioclim$bio_15 %>% terra::extract(pts_alphadiv) %>% dplyr::select(2) %>% deframe()

test_alphadiv_bioclim <- data.frame(ALPHA=rep(c("richness","shannon","Frich","Fmpd","Prich","Pmpd",
                                                "cwm_height","cwm_blade","cwm_nutlet","lifeform_index"), 5),
                                    BIOCLIM=c( rep('bio_1',10), rep('bio_7',10), rep('bio_12',10), rep('bio_15',10),
                                               rep('y',10)),
                                    F_=NA, ESS=NA, DF=NA, corr=NA, P_value=NA, B_=NA)

# transform latitude to reflect distance to the equator
alphadiv$y[alphadiv$y<0] <- alphadiv$y[alphadiv$y<0] * c(-1)

for (i in 1:nrow(test_alphadiv_bioclim)) {
  
  temp <- alphadiv[,c(test_alphadiv_bioclim$ALPHA[i], test_alphadiv_bioclim$BIOCLIM[i], 'x', 'y')] %>% na.omit()
  
  # thin to reduce computation time (as it is random the result is approx. the same)
  temp <- temp[sample(1:nrow(temp), 10000),]
  # check 'r_habit'
  if (test_alphadiv_bioclim$ALPHA[i]=="r_habit") {
    temp$r_habit[temp$r_habit=='Inf'] <- NA
    temp <- na.omit(temp)
  }

  modtest1 <- modified.ttest(temp[,test_alphadiv_bioclim$ALPHA[i]],
                             temp[,test_alphadiv_bioclim$BIOCLIM[i]],
                             temp[,c('x','y')], nclass=13)
  
  test_alphadiv_bioclim$F_[i] <- modtest1$Fstat
  test_alphadiv_bioclim$ESS[i] <- modtest1$ESS
  test_alphadiv_bioclim$DF[i] <- modtest1$dof
  test_alphadiv_bioclim$corr[i] <- modtest1$corr
  test_alphadiv_bioclim$P_value[i] <- modtest1$p.value
  test_alphadiv_bioclim$B_[i] <- lm(temp[,1] ~ temp[,2])$coefficients[2]
  
  print(paste('--- ', round(i/nrow(test_alphadiv_bioclim),2)*100,'% ---', sep=''))
  
}

# write.csv(test_alphadiv_bioclim, 'results/test_exp_alphadiv_bioclim.csv')


