

# compare observed vs expected richness at the country level to identify darkspots

source('scripts/1) import and prepare data.R')
library(ggrepel)
library(ggpubr)


# botanical countries as study unit
bot_countries <- terra::vect('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/botanicalcountries/level3.shp') %>%
  terra::project('+proj=eqearth') %>%
  terra::crop(terra::ext(c(-16921846.854212, 17080820.2067624, -6592927.59846645, 7315706.97886572)))

# plot(bot_countries)

# Magadan produces duplicated rows, it is an artifact
bot_countries <- bot_countries %>% terra::subset(bot_countries$LEVEL3_NAM!='Magadan')

# get observed richness
darkspots <- terra::extract(bot_countries, pts_observed)
darkspots$epithet <- pts_observed$epithet

# retain one observation per species and botanical country
darkspots <- darkspots %>% dplyr::select(LEVEL3_NAM, LEVEL3_COD, epithet) %>% unique() %>%
  group_by(LEVEL3_NAM, LEVEL3_COD) %>%
  summarise(observed=n()) %>%
  na.omit()



# calculate average expected probability at points of presence
richness_50km <- terra::rast('results/maps/scleria_richness_50km.tiff') # cumulative prob. x cell
exp_at_pres <- terra::extract(richness_50km, pts_observed, cells=T, ID=F, xy=T)
colnames(exp_at_pres)[colnames(exp_at_pres)=='scleria_richness'] <- 'richness'
exp_at_pres$epithet <- pts_observed$epithet # add epithet
exp_at_pres <- unique(exp_at_pres) # retain one point per species and pixel for real estimates of obs richness

# retain observed richness and expected richness at presence points
exp_at_pres <- exp_at_pres %>% group_by(cell, x, y) %>%
  summarise(observed_richness=n(),
            expected_richness=unique(richness))

# calculate residual and rasterise for downstream analyses
exp_at_pres$residual <- exp_at_pres$expected_richness-exp_at_pres$observed_richness
residual_r <- terra::rast(exp_at_pres[,c('x','y','residual')]); crs(residual_r) <- '+proj=eqearth'
# plot(residual_r)
# lines(bot_countries)



# calculate area standardized cumulative probability x botanical region, this avoids to set random thresholds
# calculate average residual per botanical country
darkspots$cum_prob <- NA
darkspots$areakm <- NA
darkspots$residual_mean <- NA
darkspots$residual_sd <- NA
for (i in 1:nrow(darkspots)) {
  # select polygon
  aoi_l <- bot_countries[bot_countries$LEVEL3_COD==darkspots$LEVEL3_COD[i]]
  darkspots$areakm[i] <- terra::expanse(aoi_l, unit='km')
  # stand prob
  aoi_r <- richness_50km %>% terra::crop(aoi_l, mask=T)
  darkspots$cum_prob[i] <- aoi_r %>% as.data.frame() %>% deframe() %>% sum()
  # residual
  aoi_r <- residual_r %>% terra::crop(aoi_l, mask=T)
  darkspots$residual_mean[i] <- aoi_r %>% as.data.frame() %>% deframe() %>% mean(na.rm=T)
  darkspots$residual_sd[i] <- aoi_r %>% as.data.frame() %>% deframe() %>% sd(na.rm=T)
}

# observed species per km2
darkspots$obs_st_1000 <- darkspots$observed / darkspots$areakm * 1000

# expected richness per km2
darkspots$exp_st_1000 <- darkspots$cum_prob /darkspots$areakm * 1000



# plot
# head(darkspots)
model <- lm(residual_mean ~ observed, darkspots[,c('residual_mean','observed')])
anova(lm(residual_mean~observed, darkspots))
pred_df <- predict(model, newdata=darkspots[,c('residual_mean','observed')], interval="confidence", level=0.95)
darkspots$fit <- pred_df[,"fit"]
darkspots$lwr <- pred_df[,"lwr"]
darkspots$upr <- pred_df[,"upr"]
# calculate distance to estimate
darkspots$deviation <- darkspots$residual_mean - darkspots$fit
# calculate distance to estimate
darkspots$upr_true <- darkspots$residual_mean > darkspots$fit


# identify top 10 darkspots for label
darkspots$LEVEL3_COD_2 <- darkspots$LEVEL3_COD
table(darkspots$deviation<2.9)
darkspots$LEVEL3_COD_2[which(darkspots$deviation<2.9)] <- NA


# retrieve information of level1 areas
wgsrpd <- merge(terra::vect('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/botanicalcountries/level1.shp') %>% as.data.frame(), terra::vect('C:/Users/javie/OneDrive/ACADEMICO/proyectos/world_rasters/botanicalcountries/level3.shp') %>% as.data.frame()) %>% dplyr::select(LEVEL1_NAM, LEVEL3_NAM, LEVEL3_COD)
darkspots <- left_join(darkspots, wgsrpd, by=c('LEVEL3_NAM','LEVEL3_COD'))
darkspots <- darkspots[darkspots$LEVEL1_NAM!='PACIFIC',]


g1 <- ggplot(aes(y=residual_mean, x=observed, size=areakm, color=residual_mean), data=darkspots) +
  geom_point() +
  xlab('Observed richness') + ylab(NULL) +
  geom_text_repel(aes(label=LEVEL3_COD_2), size=3, nudge_x=2, nudge_y=0.5, col='blue') +
  geom_smooth(method='lm', level=0.95, col='red', se=T) +
  geom_abline(slope=0, intercept=0, linetype='dashed') +
  scale_y_continuous(limits=c(-2.5,10)) + scale_x_continuous(limits=c(2,45)) +
  theme_bw() +
  theme(legend.position ='right', legend.title=element_text(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=9, color='black')) +
  labs(colour="Deviation", size='Area (km2)')


# darkspots$LEVEL1_NAM <- factor(darkspots$LEVEL1_NAM, levels=c("SOUTHERN AMERICA","AFRICA","ASIA-TROPICAL","AUSTRALASIA","NORTHERN AMERICA","ASIA-TEMPERATE"))
# levels(darkspots$LEVEL1_NAM) <- c("SAM","AFR","ASTR","AUS","NAM","ASTE")

g2 <- ggplot(aes(y=residual_mean, x=LEVEL1_NAM), data=darkspots) +
  geom_boxplot() +
  xlab(NULL) + ylab('Expected - observed richness') +
  theme_bw() +
  geom_abline(slope=0, intercept=0, linetype='dashed') +
  scale_y_continuous(limits=c(-2.5,10)) +
  theme(legend.position ='right', legend.title=element_text(),
        axis.text.x=element_text(size=10, color='black'),
        plot.margin = margin(t=5, r=5, b=18, l=20)) +
  labs(colour="Deviation", size='Area (km2)')


ggarrange(g2, g1, widths=c(1, 2.2))


