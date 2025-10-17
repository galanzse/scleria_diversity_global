

# import and explore trait data

scl_traits <- data_final[['traits']] %>%
  merge(scl_taxonomy[,c('scientific_name','epithet')], by='scientific_name') %>%
  dplyr::select(-scientific_name)
scl_traits <- merge(scl_traits, scl_taxonomy[,c('epithet','section')])

# impute by section x life form means
for (i in 1:nrow(scl_traits)) {
  for (t in c('life_form_simp', 'height', 'blade_area', 'nutlet_volume')) {
    if (is.na(scl_traits[i,t])) {
      temp <- scl_traits %>%
        subset(section==scl_traits$section[i] & life_form_simp==scl_traits$life_form_simp[i]) %>%
        dplyr::select(all_of(t)) %>% colMeans(na.rm=T)
      scl_traits[i,t] <- temp
    }
  }
} 

scl_traits[scl_traits=='NaN'] <- NA # correct format

scl_traits[scl_traits$epithet=='khasiana',c('height', 'blade_area', 'nutlet_volume')] <- scl_traits %>% subset(section=='Elatae') %>% dplyr::select(height, blade_area, nutlet_volume) %>% colMeans(na.rm=T) # S. khasiana is imputed from section data only

# transform and scale (Mammola et al. 2021)
scl_traits$height <- scl_traits$height %>% log() %>% scale(center=F) %>% as.vector()
scl_traits$blade_area <- scl_traits$blade_area %>% log() %>% scale(center=F) %>% as.vector()
scl_traits$nutlet_volume <- scl_traits$nutlet_volume %>% log() %>% scale(center=F) %>% as.vector()

# hist(scl_traits$height)
# hist(scl_traits$blade_area)
# hist(scl_traits$nutlet_volume)

# prepare matrix for analyses
scl_traits <- scl_traits %>% dplyr::select(epithet, life_form_simp, height, blade_area, nutlet_volume)

# species to work with (present in all datasets and with infrageneric information)
v_species <- intersect(intersect(scl_occurrences$epithet, scl_traits$epithet), scl_taxonomy$epithet)
scl_occurrences <- scl_occurrences %>% filter(epithet %in% v_species) # filter occurrences without trait data
pts_observed <- scl_occurrences %>% vect(geom=c('x','y'), 'epsg:4326') %>% project('+proj=eqearth') # fix
scl_traits <- scl_traits %>% filter(epithet %in% v_species) # fix
scl_traits$life_form_simp <- as.numeric(as.factor(scl_traits$life_form_simp)) # fix

scl_taxonomy <- scl_taxonomy %>% subset(epithet%in%v_species)

scl_traits <- scl_traits %>% dplyr::select(epithet, life_form_simp, height, blade_area, nutlet_volume) 
rownames(scl_traits) <- scl_traits$epithet; scl_traits$epithet <- NULL

# We will use a dendrogram so there are values of Frich for all species
str(scl_traits)
scl_dist <- cluster::daisy(as.matrix(scl_traits), metric='gower', type=list('factor'=1,'numeric'=2:4))
scl_dendrogram <- hclust(scl_dist, method='average') # clustering (UPGMA)

# # plot
# scl_dend_phy <- ape::as.phylo(scl_dendrogram)
# 
# # Pick colors for each group
# group_colors <- c('Scleria'='forestgreen', 'Hypoporum'='orange', 'Trachylomia'='violet', 'Browniae'='blue')
# 
# # Map species to colors
# tip_colors <- group_colors[scl_taxonomy$subgenus]
# 
# # Make sure the order matches tree$tip.label
# names(tip_colors) <- scl_taxonomy$epithet
# tip_colors <- tip_colors[scl_dend_phy$tip.label]
# plot.phylo(scl_dend_phy, type='fan', tip.color=tip_colors, cex=0.75, main='Functional dendrogram')
# legend("center", y.intersp=0.3, x.intersp=0.5,
#        legend=names(group_colors), col=group_colors,
#        pch=19, pt.cex=1.2, bty="n",  title="")



