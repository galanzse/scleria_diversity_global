

# OBSERVED ALPHA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY, AND SES


source('scripts/1) import and prepare data.R')

library(BAT) #Frich


# DATA
Xnull=499
scl_traits # traits, distance and dendrogram
scl_dist
scl_dendrogram
imputed_trees # imputed phylogenies to account for uncertainty
observed_assemblages # community matrix



# NULL COMMUNITIES ####

# we are going to define one global pool: any species can occur in any grid (Swenson et al. 2012)
# null.model=c('richness','independentswap'). Makes sense to compare to a similar richness structure

# number of null communities computed
Nswaps=10000

# # observed assemblages
# observed_assemblages[1:5, 1:5]
# # randomize matrix and save
# null_obs_ass <- list()
# for (i in 1:Xnull) {
#   # randomize assemblage
#   null_obs_ass[[i]] <- picante::randomizeMatrix(samp=observed_assemblages[,-c(1:3)],
#                                                 null.model='richness') # iterations=Nswaps
#    # progress
#    print(paste(round(i/Xnull, 4)*100, '%'))
# }
# save(null_obs_ass, file="results/null_ass/null_obs_ass_rich.RData")



# CALCULATE ALPHA DIVERSITY ####

# calculate Faith and Pmpd as the mean values across N imputed trees to account for imputation bias
Ntrees <- 10

# site x species presence matrix to subset communities from
observed_assemblages

# results df
obs_indices <- observed_assemblages[,c('cell','x','y')]


# richness
obs_indices$richness <- rowSums(observed_assemblages[,-c(1:3)])


# functional composition
temp <- BAT::cwm(observed_assemblages[,-c(1:3)],
                 scl_traits[,c('height','blade_area','nutlet_volume')],
                 abund=FALSE, na.rm=FALSE)
obs_indices$cwm_height <- temp[,'height']
obs_indices$cwm_blade <- temp[,'blade_area']
obs_indices$cwm_nutlet <- temp[,'nutlet_volume']


# life form: proportion of annual species
obs_indices$n_ann <- NA
v_annuals <- scl_traits %>% subset(life_form_simp==1) %>% rownames()
for (i in 1:nrow(obs_indices)) {
  # select assemblage
  temp <- observed_assemblages[i,-c(1:3)]
  # retain present taxa
  temp <- names(temp[,colSums(temp)==1])
  # get proportion of annual species
  obs_indices$n_ann[i] <- length(which(temp %in% v_annuals))
  # progress
  print(i)
}

# functional indexes
obs_indices$Frich <- BAT::alpha(observed_assemblages[,-c(1:3)], scl_dendrogram)[,'Richness']
obs_indices$Fmpd <- picante::mpd(observed_assemblages[,-c(1:3)],
                                 as.matrix(cophenetic(scl_dendrogram)), abundance.weighted=FALSE)


# phylogenetic diversity
# mat_Prich <- matrix(nrow=nrow(observed_assemblages), ncol=Ntrees)
# mat_Pmpd <- mat_Prich
# for (i in 1:Ntrees) {
#   mat_Prich[,i] <- BAT::alpha(observed_assemblages[,-c(1:3)], sample(imputed_trees,1)[[1]])
#   mat_Pmpd[,i] <- picante::mpd(observed_assemblages[,-c(1:3)],
#                                as.matrix(cophenetic(sample(imputed_trees,1)[[1]])),
#                                abundance.weighted=FALSE)
#   print(i)
# }
# 
# # save matrices
# obs_mat_Pdiv <- list()
# obs_mat_Pdiv[['Prich']] <- mat_Prich
# obs_mat_Pdiv[['Pmpd']] <- mat_Pmpd
# save(obs_mat_Pdiv, file="results/obs_mat_Pdiv.RData")
load("results/obs_mat_Pdiv.RData")

# calculate mean and save
obs_indices$Prich <- rowMeans(obs_mat_Pdiv$Prich, na.rm=F)
obs_indices$Pmpd <- rowMeans(obs_mat_Pdiv$Pmpd, na.rm=F)


# checkpoint
# write.table(obs_indices, 'results/obs_alpha_div.txt')
obs_indices <- read.csv("results/obs_alpha_div.txt", sep="")

# pairs(obs_indices[,c("richness","Frich","Fmpd","Prich","Pmpd")],
#       upper.panel=NULL, main='pairs(observed alpha diversity)')



# NULL MODELS: compute indices across null communities to test for significancy ####

# # random communities via independentswap
# load("results/null_obs_ass.RData")
# head(null_obs_ass[[1]])
# 
# alpha_obs_null <- matrix(nrow=nrow(observed_assemblages), ncol=Xnull) # results for every combination
# rownames(alpha_obs_null) <- rownames(observed_assemblages)
# 
# null.Frich <- alpha_obs_null
# null.Fmpd <- alpha_obs_null
# null.Prich <- alpha_obs_null
# null.Pmpd <- alpha_obs_null
# 
# # check match
# length(scl_dendrogram$labels) == length(colnames(null_obs_ass[[1]]))
# table(scl_dendrogram$labels %in% colnames(null_obs_ass[[1]]))
# table(rownames(obs_indices)==rownames(null_obs_ass[[1]]))
# 
# for (i in 1:Xnull) {
# 
#   # Frich
#   null.Frich[,i] <- BAT::alpha(null_obs_ass[[i]], scl_dendrogram)
# 
#   # Fmpd
#   null.Fmpd[,i] <- picante::mpd(null_obs_ass[[i]], as.matrix(cophenetic(scl_dendrogram)), abundance.weighted=FALSE)
# 
#    # temporal matrices to store phylo diversity
#    mat_Prich <- matrix(nrow=nrow(null_obs_ass[[i]]), ncol=Ntrees)
#    mat_Pmpd <- matrix(nrow=nrow(null_obs_ass[[i]]), ncol=Ntrees)
# 
#   for (t in 1:Ntrees) {
#     mat_Prich[,t] <- BAT::alpha(null_obs_ass[[i]], sample(imputed_trees,1)[[1]])
#     mat_Pmpd[,t] <- picante::mpd(null_obs_ass[[i]],
#                                  as.matrix(cophenetic(sample(imputed_trees,1)[[1]])),
#                                  abundance.weighted=FALSE)
#     print(paste(t))
#   }
# 
#   null.Prich[,i] <- rowMeans(mat_Prich, na.rm=TRUE)
#   null.Pmpd[,i] <- rowMeans(mat_Pmpd, na.rm=TRUE)
# 
# 
#   # progress
#   print(paste('--- ', round(i/Xnull*100, 2), '% ---', sep=''))
# 
# }
# 
# # store in list
# alpha_obs_null_list <- list() # save
# alpha_obs_null_list[['Frich']] <- null.Frich
# alpha_obs_null_list[['Fmpd']] <- null.Fmpd
# alpha_obs_null_list[['Prich']] <- null.Prich
# alpha_obs_null_list[['Pmpd']] <- null.Pmpd
# 
# # save
# save(alpha_obs_null_list, file='results/alpha_obs_null_list.Rdata')
load("results/alpha_obs_null_list.RData")



# identify significant richness and diversity patterns
obs_indices$perc_Frich <- NA
obs_indices$perc_Prich <- NA
obs_indices$perc_Fmpd <- NA
obs_indices$perc_Pmpd <- NA

for (i in 1:nrow(obs_indices)) {
  if (obs_indices$richness[i]>1) {
    obs_indices$perc_Prich[i] <- ecdf(alpha_obs_null_list$Prich[i,])(obs_indices$Prich[i])
    obs_indices$perc_Frich[i] <- ecdf(alpha_obs_null_list$Frich[i,])(obs_indices$Frich[i])
    obs_indices$perc_Fmpd[i] <- ecdf(alpha_obs_null_list$Fmpd[i,])(obs_indices$Fmpd[i])
    obs_indices$perc_Pmpd[i] <- ecdf(alpha_obs_null_list$Pmpd[i,])(obs_indices$Pmpd[i])
  }
  # progress
  print(i)
}

# table(obs_indices$perc_Prich>0.95)/length(obs_indices$perc_Prich)

# # save
# write.table(obs_indices, 'results/obs_alpha_div.txt')


