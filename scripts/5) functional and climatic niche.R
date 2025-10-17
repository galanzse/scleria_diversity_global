

# EXPLORE RELATIONSHIP BETWEEN FUNCTIONALITY AND CLIMATIC NICHE IN SCLERIA AT THE SUBGENUS LEVEL

source('scripts/1) import and prepare data.R')


library(geiger)
library(terra)

# library(randomForest)
# library(ggpubr)
# library(caret)
# library(caTools)
# library(SpatialPack) # modified.ttest

load("results/imputed_trees.RData")



# FUNCTIONAL DIFFERENCES ####

# use original data to avoid imputation
LHS_data <- merge(data_final[['traits']], scl_taxonomy, all.x=T) %>%
  dplyr::select('epithet',"subgenus","section","life_form","life_form_simp","height","blade_area","nutlet_volume")

# format data to run phylogenetic anova
temp <- LHS_data[,c('epithet','subgenus',"nutlet_volume")] %>%
  # avoid imputation
  na.omit() %>%
  # exclude species without occurrence data
  subset(epithet %in%scl_occurrences$epithet)

dat <- temp$nutlet_volume
names(dat) <- temp$epithet
dat <- log(dat)

group <- temp$subgenus
names(group) <- temp$epithet
group <- as.factor(group)


geiger::aov.phylo(dat ~ group, phy=imputed_trees[[1]], nsim=1000)
# phytools::phylANOVA(tree=keep.tip(imputed_trees[[1]], temp$epithet), x=temp$subgenus, y=temp$nutlet_volume, posthoc=TRUE)


# 
chisq.test(table(LHS_data$life_form_simp, LHS_data$subgenus))

table(LHS_data$life_form_simp, LHS_data$subgenus)
colSums(table(LHS_data$life_form_simp, LHS_data$subgenus))


# plot
v_traits <- c('height','blade_area','nutlet_volume')
LHS_data[,v_traits] <- LHS_data[,v_traits] %>% apply(2, log)

long_LHS_data <- LHS_data %>%
  dplyr::select(subgenus, all_of(v_traits)) %>%
  pivot_longer(2:4, names_to='trait', values_to='value') %>%
  na.omit()

long_LHS_data$subgenus <- as.factor(long_LHS_data$subgenus)
long_LHS_data$subgenus <- factor(long_LHS_data$subgenus, levels=c("Hypoporum", "Trachylomia", "Browniae", "Scleria"))

# there is a clear functional segregation at the subgenus level for leaf and nutlet traits
long_LHS_data$trait <- factor(long_LHS_data$trait, c('height','blade_area','nutlet_volume'))
lbs = setNames(c("'log Maximum height (cm)'",
                 "'log Blade area ('*cm^2*')'", 
                 "'log Nutlet volume ('*mm^3*')'"),
               c('height','blade_area','nutlet_volume'))

ggplot(aes(x=subgenus, y=value), data=long_LHS_data) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) +
  facet_wrap(.~trait, scales="free", ncol=3, labeller=as_labeller(lbs, label_parsed)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=-45, vjust=1.2, hjust=0))



# CLIMATIC DIFFERENCES ####

scl_bioclim <- terra::rast('SDM/world_rasters_rs2.5/worldclim.tiff')[[c("bio_1","bio_7","bio_12","bio_15")]]
scl_bioclim <- scl_bioclim %>% terra::extract(scl_occurrences[,c('x','y')], ID=F)
scl_bioclim$epithet <- scl_occurrences$epithet
scl_bioclim <- merge(scl_bioclim, scl_taxonomy[,c('epithet','section','subgenus')], all.x=T)

scl_bioclim <- scl_bioclim %>% dplyr::select(epithet, subgenus, bio_1, bio_7, bio_12, bio_15) %>%
  group_by(epithet, subgenus) %>%
  summarise(bio_1=mean(bio_1, na.rm=T),
            bio_7=mean(bio_7, na.rm=T),
            bio_12=mean(bio_12, na.rm=T),
            bio_15=mean(bio_15, na.rm=T))

# this line avoids the error: Error en pic(x, td$phy): 'phy' is not rooted and fully dichotomous (?)
scl_bioclim <- scl_bioclim %>% subset(epithet %in% temp$epithet)

# no need to log-trans for anovas
dat <- scl_bioclim$bio_15
names(dat) <- scl_bioclim$epithet

group <- scl_bioclim$subgenus
names(group) <- scl_bioclim$epithet
group <- as.factor(group)

geiger::aov.phylo(dat ~ group, phy=imputed_trees[[1]], nsim=1000)


# plot
long_scl_bioclim <- scl_bioclim %>% pivot_longer(3:6, names_to='bioclim', values_to='value') %>% na.omit()
long_scl_bioclim$bioclim <- as.factor(long_scl_bioclim$bioclim)
levels(long_scl_bioclim$bioclim) <- c('Annual Mean Temperature', 'Annual Precipitation', 'Precipitation Seasonality','Temperature Annual Range')
long_scl_bioclim$bioclim <- factor(long_scl_bioclim$bioclim, c('Annual Mean Temperature', 'Temperature Annual Range', 'Annual Precipitation', 'Precipitation Seasonality'))


long_scl_bioclim$subgenus <- as.factor(long_scl_bioclim$subgenus)
long_scl_bioclim$subgenus <- factor(long_scl_bioclim$subgenus, levels=c("Hypoporum", "Trachylomia", "Browniae", "Scleria"))


lbs = setNames(c("'Annual Mean Temperature ' (degree*C)",
                 "'Temperature Annual Range ' (degree*C)",
                 "'Annual Precipitation (mm)'",
                 "'Precipitation Seasonality (CV)'"),
               levels(long_scl_bioclim$bioclim))

ggplot(aes(x=subgenus, y=value), data=long_scl_bioclim) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) +
  facet_wrap(.~bioclim, scales="free", nrow=1, labeller=as_labeller(lbs, label_parsed)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=-45, vjust=1.2, hjust=0))



# DENDROGRAM COMPARISON ####

# climatic dendrogram
temp <- scl_bioclim
temp$bio_1 <- temp$bio_1 %>% log() %>% scale(center=F) %>% as.vector()
temp$bio_7 <- temp$bio_7 %>% log() %>% scale(center=F) %>% as.vector()
temp$bio_12 <- temp$bio_12 %>% log() %>% scale(center=F) %>% as.vector()
temp$bio_15 <- temp$bio_15 %>% log() %>% scale(center=F) %>% as.vector()

temp <- as.matrix(temp[,c('bio_1','bio_7','bio_12','bio_15')])
rownames(temp) <- scl_bioclim$epithet

scl_dist <- cluster::daisy(as.matrix(temp), metric="gower")
clim_dendrogram <- hclust(scl_dist, method="average") # clustering (UPGMA)


# trait dendrogram
scl_dist <- cluster::daisy(as.matrix(scl_traits[clim_dendrogram$labels,]),
                           metric="gower", type=list('factor'=1,'numeric'=2:4))
trait_dendrogram <- hclust(scl_dist, method="average") # clustering (UPGMA)


# install.packages('dendextend')
# library(dendextend)
# dendextend::entanglement(trait_dendrogram, clim_dendrogram)
# dendextend::cor_cophenetic(trait_dendrogram, clim_dendrogram)


