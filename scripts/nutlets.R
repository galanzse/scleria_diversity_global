

# IS NUTLET SIZE CORRELATED TO DISPERSAL CAPACITY AT THE SPECIES LEVEL?


source('scripts/1) import and prepare data.R')

library(readxl)
library(geiger)
library(ggpubr)


# nutlets
scl_nutlets <- read_excel("data/nutlets.xlsx")
scl_nutlets$epithet <- sapply(strsplit(scl_nutlets$scientific_name, " "), function(x) x[2])

# retain nutlets in tree for phylogenetic ANOVAs
scl_nutlets <- scl_nutlets %>% subset(epithet %in% imputed_trees[[1]]$tip.label) %>% dplyr::select(-length, -width, -scientific_name)
scl_nutlets <- merge(scl_nutlets, scl_taxonomy[,c('epithet','subgenus','section')])

# ranges
scl_ranges <- terra::rast('results/maps/scleria_probocc_50.tiff')
names(scl_ranges) <- sapply(strsplit(names(scl_ranges), " "), function(x) x[2])
names(scl_ranges)[names(scl_ranges)=='baroni-clarkei'] <- 'baroniclarkei'
names(scl_ranges)[names(scl_ranges)=='flagellum-nigrorum'] <- 'flagellumnigrorum'
names(scl_ranges)[names(scl_ranges)=='novae-hollandiae'] <- 'novaehollandiae'

# tree
load("results/imputed_trees.RData")


# distribution range size
scl_nutlets$n_obs <- NA
scl_nutlets$EOO <- NA

# threshold to convert probability of occurrence to presence
thres1=0.8

for (i in 1:nrow(scl_nutlets)) {
  
  # presence data
  occ_obs <- scl_occurrences[scl_occurrences$epithet==scl_nutlets$epithet[i],] %>% vect(geom=c('x','y'), 'epsg:4326')
  
  # expected localities with very high probability
  occ_exp <- scl_ranges[[scl_nutlets$epithet[i]]] %>% as.data.frame(xy=T)
  colnames(occ_exp) <- c('x','y','prob')
  occ_exp <- occ_exp[occ_exp$prob>thres1,]
  
  if (nrow(occ_exp)==0) {
    scl_nutlets$n_obs[i] <- nrow(occ_obs)
    scl_nutlets$EOO[i] <- occ_obs %>% terra::convHull() %>% terra::expanse(unit='km')
  } else {
    occ_exp <- occ_exp %>% vect(geom=c('x','y'), '+proj=eqearth') %>%
      terra::project('epsg:4326') %>% rbind(occ_obs)
    scl_nutlets$n_obs[i] <- nrow(occ_exp)
    scl_nutlets$EOO[i] <- occ_exp %>% terra::convHull() %>% terra::expanse(unit='km')
  }
  
  print(i)
}


# categorize size
qnt1 <- quantile(scl_nutlets$volume, c(0,1/3,2/3,1), na.rm=T)[2:3]
scl_nutlets$volume_cat <- cut(scl_nutlets$volume, breaks=c(-Inf, qnt1[1], qnt1[2], Inf), labels=c("Small","Medium","Big"))

# categorize shape
qnt1 <- quantile(scl_nutlets$shape, c(0,1/3,2/3,1), na.rm=T)[2:3]
scl_nutlets$shape_cat <- cut(scl_nutlets$shape, breaks=c(-Inf, qnt1[1], qnt1[2], Inf), labels=c("Round","Intermediate","Elongated"))


# plots
boxplot(log(scl_nutlets$EOO+1) ~ scl_nutlets$type_ornamentation)
boxplot(log(scl_nutlets$EOO+1) ~ scl_nutlets$ornamented)
boxplot(log(scl_nutlets$EOO+1) ~ scl_nutlets$hairy)
boxplot(log(scl_nutlets$EOO+1) ~ scl_nutlets$volume_cat)

plot(log(scl_nutlets$EOO+1) ~ scl_nutlets$volume)
plot(log(scl_nutlets$EOO+1) ~ scl_nutlets$weight)
plot(log(scl_nutlets$EOO+1) ~ scl_nutlets$shape)

plot(log(scl_nutlets$weight) ~ scl_nutlets$volume)


# importance of evolutionary relationships in EOO
boxplot(log(scl_nutlets$EOO+1) ~ scl_nutlets$subgenus)
boxplot(log(scl_nutlets$EOO+1) ~ scl_nutlets$section, las=2, xlab=NULL, cex.axis=0.7)


# phylogenetic anova
temp <- scl_nutlets %>% dplyr::select(epithet, EOO, volume_cat) %>% na.omit()
dat <- log(temp$EOO+1)
names(dat) <- temp$epithet
group <- temp$volume_cat
names(group) <- temp$epithet
geiger::aov.phylo(dat ~ group, keep.tip(imputed_trees[[1]], temp$epithet))


# my_comparisons <- list( c("Small", "Medium"), c("Small", "Big"), c("Medium", "Big"))
ggplot(aes(y=log(EOO), x=volume), data=scl_nutlets) + # AOO
  geom_point() +
  ylab('log(EOO)') + xlab('Nutlet size') +
  theme_classic() 
# stat_compare_means(method='t.test', paired=F, comparisons=my_comparisons, bracket.size=.1,
#                    label="p.signif", hide.ns=T, vjust=0.4)

ggplot(aes(y=log(EOO), x=subgenus), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('log(EOO)') + xlab('') +
  theme_classic() 



anova(lm(log(EOO) ~ nutlet_cat, data=scl_nutlets))
ggplot(aes(x=nutlet_cat, y=log(EOO)), data=scl_nutlets) + # AOO
  geom_boxplot() +
  geom_smooth(method='lm') +
  xlab('height') + ylab('log(EOO)') +
  theme_classic()

ggplot(aes(y=y_range, x=nutlet_cat), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('Latitudinal range') + xlab('Nutlet size category') +
  theme_classic()

ggplot(aes(y=x_range, x=nutlet_cat), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('Longitudinal range') + xlab('Nutlet size category') +
  theme_classic()

