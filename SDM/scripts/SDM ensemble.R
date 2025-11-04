

# Run ensemble SDMs after Gonzalez-Moreno (Xylella) for species with >= 10 observations
source('SDM/scripts/SDM organize analyses.R')
source('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/functions/fun_PAtable.R')


library(biomod2)
# library(tidyterra)


# species to model
# scl_occ_ens <- scl_sdm_alg %>% subset(algorithm == 'ensemble')


# parameters
n_pa_sets = 3 # number of bg sets to reduce bias
dtp = 500000 # distance to presences to define target_raster

# model selection
v_models <- c('MAXENT','RF','GBM','ANN','GLM','GAM')

# if path is too long the function yields error, plus there are too many output files, so lets work in the desktop
setwd('C:/Users/javie/Desktop/biomod_output')

# plot parameters
cuts <- seq(0, 1, 0.1)
pal <- colorRampPalette(c("white","black"))


# library(readxl)
scl_occ_ens <- read_excel("C:/Users/javie/Desktop/redo_ens.xlsx")

for (i in 24:nrow(scl_occ_ens)) {
  
  # S. polycarpa yields error so it needs to be done separately (already done)
  if (scl_occ_ens$species[i]=='Scleria polycarpa Boeckeler') next
  
  
  # DATA PREPARATION: importar una especie
  df_p <- data_final$occurrences %>% subset(scientific_name==scl_occ_ens$species[i])
  pts_p <- df_p %>% terra::vect(geom=c('x','y'), crs='+proj=longlat')
  
  # model in ecoregions that intersect the species EOO
  aoi_eco <- AOI %>% terra::mask(terra::convHull(pts_p)) %>% terra::aggregate(dissolve=TRUE)
  
  # crop predictors
  r <- crop_unc_rst_predictors_2.5min %>% terra::crop(aoi_eco, mask=T)
  
  
  # bias file: model in ecoregions were Scleria spp. has been found
  # target-group background (bg) (https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13442)
  target_points <- scl_occurrences %>% unique() %>%
    terra::vect(geom=c('x','y'), 'epsg:4326') %>%
    terra::intersect(aoi_eco)
  
  # some species have large distribution ranges. We need to restrict the generation of bg points to sites not too far to presences (e.g. 500km) so bg points are sites environmentally kind of close to the presence points
  bff_p <- pts_p %>% terra::buffer(width=dtp) %>% terra::aggregate(dissolve=TRUE) # 500km
  target_points <- target_points %>% terra::crop(bff_p)
  
  # 2d kernel density estimation based on Scleria spp. observations 500km from the species presence
  target_raster <- ks::kde(geom(target_points)[,c('x','y')]) %>%
    raster::raster() %>% terra::rast() %>%
    terra::crop(aoi_eco, mask=T) %>%
    terra::resample(terra::rast(r), method='bilinear') %>%
    terra::crop(bff_p, mask=T)

  # normalize
  target_raster <- target_raster - min(target_raster[], na.rm=T) 
  target_raster <- spatialEco::raster.transformation(target_raster, trans="norm")

  # generate df of background points
  bg <- as.data.frame(target_raster, xy=T) %>% subset(layer>0.1)
  colnames(bg)[colnames(bg)=='layer'] <- 'prob'
  
  # plot(target_raster)
  # points(bg[,c('x','y')], col='orange')
  # points(pts_p, col='red')
  
  
  # prepare pseudoabsence data to include personal selection of bg from bias file
  df_pres <- df_p[,c('x','y')]
  df_pres$origin <- 1 # 1 indicates true presence
  
  # if there are n_pa*n_pa_sets bg points, then define n_pa_sets independent sets of n_pa observations
  n_pa = 1000
  if (nrow(bg)>=c(n_pa*n_pa_sets)) {
    n_pa = 1000 # PA size
    df_pseab <- bg[sample(1:nrow(bg), n_pa*n_pa_sets, replace=F),c('x','y')]
    df_pseab$origin <- 0 # indicates pseudoabsence
    df_respvar <- rbind(df_pres, df_pseab) # merge presences with PAsets
    pts_respvar <- terra::vect(df_respvar, geom=c('x','y'), crs='epsg:4326')
  # if there are less than n_pa*n_pa_sets bg points, just split the bg points in three sets
  } else {
    n_pa <- trunc(nrow(bg)/n_pa_sets) # PA size
    df_pseab <- bg[1:c(n_pa*n_pa_sets),c('x','y')] # make sure nrow(df_pseab) == n_pa*n_pa_sets
    df_pseab <- df_pseab[sample(1:nrow(df_pseab)),] # randomize to avoid geographic clustering
    df_pseab$origin <- 0 # indicates pseudoabsence
    df_respvar <- rbind(df_pres, df_pseab) # merge presences with PAsets
    pts_respvar <- terra::vect(df_respvar, geom=c('x','y'), crs='epsg:4326')
  }
  
  # generate PA.user.table to interpret input of resp.var: this table allows the function to interpret pts_respvar
  df_sets <- PAtable(npres=nrow(df_p), PA.nb.rep=n_pa_sets, PA.nb.absences=n_pa)
  
  
  # BUILD MODELS
  
  # format data
  myBiomodData <- biomod2::BIOMOD_FormatingData(
    resp.name = scl_occ_ens$species_simp[i],
    resp.var = pts_respvar,
    expl.var = r,
    filter.raster = TRUE,
    dir.name = getwd(),
    PA.nb.rep = n_pa_sets,
    PA.strategy = 'user.defined',
    PA.user.table = df_sets,
    eval.resp.var = NULL, # presencias para validacion/ necesitaria ausencias reales
    na.rm = TRUE)
  
  
  # build individual models
  myBiomodModelOut <- biomod2::BIOMOD_Modeling(
    bm.format=myBiomodData,
    models = v_models, #  https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.13555 
    OPT.strategy = "bigboss", # defined by biomod2 team
    CV.strategy = 'random',
    CV.nb.rep = 5, # 5 kfold cross validation
    CV.perc = 0.8, # https://github.com/biomodhub/biomod2/issues/546
    prevalence = 0.5, # weighted sum of presences equals weighted sum of absences
    var.import = 3,
    metric.eval = c('TSS','ROC','KAPPA',"BOYCE"), # https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13515
    do.progress = TRUE)
  
  # get ALL evaluations and save (models used in ensemble are those with BOYCE>0.7)
  myBiomodModelEval <- biomod2::get_evaluations(myBiomodModelOut)
  write.table(myBiomodModelEval, file=paste('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/ensemble/evaluations/simple/', scl_occ_ens$species_simp[i], '.txt', sep=''))
  
  # get variable importance for simple models used in ensemble BOYCE>0.7 and save
  myBiomodModelImp <- biomod2::get_variables_importance(myBiomodModelOut) %>%
    subset(full.name %in% subset(myBiomodModelEval, metric.eval=='BOYCE' & validation>=0.7)$full.name)
  write.table(myBiomodModelImp, file=paste('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/ensemble/var_importance/', scl_occ_ens$species_simp[i], '.txt', sep='')) 

  
  # ensemble individual models
  myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen = "all", em.by = 'all', em.algo = 'EMwmean',
    metric.select = 'BOYCE', metric.select.thresh = 0.7,
    metric.eval = c('TSS','ROC','KAPPA',"BOYCE"),
    EMci.alpha = 0.05,
    EMwmean.decay = 'proportional')
  
  # ensemble evaluation
  myBiomodModelEmEval <- biomod2::get_evaluations(myBiomodEM) %>%
    dplyr::select(full.name, algo, metric.eval, calibration, specificity, sensitivity) %>%
    mutate(
      sensitivity = sensitivity / 100,  # percentage to proportion
      specificity = specificity / 100,  # percentage to proportion
      Youdens_Index = sensitivity + specificity - 1) %>% # Youden's I
    arrange(desc(Youdens_Index))
  
  # get ensemble evaluation and save (models used in ensemble are those with BOYCE>0.7)
  write.table(myBiomodModelEmEval, file=paste('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/ensemble/evaluations/ensemble/', scl_occ_ens$species_simp[i], '.txt', sep=''))
  

  # project model
  tryCatch({
    myBiomodEMProj <- biomod2::BIOMOD_EnsembleForecasting(
      bm.em = myBiomodEM,
      bm.proj = NULL, # we are directly projecting the ensembles
      models.chosen = 'all',
      proj.name = scl_occ_ens$species_simp[i],
      new.env = r,
      buid.clamping.mask=FALSE,
      clamp = TRUE)
  }, error=function(e){})
  

  # get raster and normalize because values are stored without decimals (0 to 1000 instead of 0 to 1)
  proj_temp <- myBiomodEMProj %>% get_predictions()
  names(proj_temp) <- scl_occ_ens$species_simp[i]
  proj_temp <- proj_temp/1000
  
  # save
  writeRaster(proj_temp, paste('C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/predictions/', scl_occ_ens$species[i], '.tiff', sep=''), overwrite=T)
  
  
  # create map
  path1 <- 'C:/Users/javie/OneDrive/ACADEMICO/proyectos/scleria/scleria global diversity/SDM/results/maps/'
  png(paste(path1, scl_occ_ens$species[i], '.png', sep=''), bg="white")
  plot(proj_temp, breaks=cuts, col=pal(10), main=scl_occ_ens$species[i])
  lines(aoi_eco, col='black')
  points(pts_p, col='green', pch=16)
  dev.off()
  
  
  # progress
  print(paste(scl_occ_ens$species[i], '--- Done!', '[', round(i/nrow(scl_occ_ens),2)*100, '% ]'))

}


# # plot variable importance
# ggplot(aes(x=algo, y=var.imp, colour=expl.var), data=myBiomodModelImp) +
#   geom_boxplot() +
#   theme_bw()


