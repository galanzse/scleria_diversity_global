


# EXPECTED (SDMs) ALPHA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY, AND SES
# RESULTS ARE ABUNDANCE WEIGHTED!


source('scripts/1) import and prepare data.R')

library(BAT) # Frich
library(vegan) # diversity(x, index="shannon")


# DATA
Xnull # 499
scl_traits # traits, distance and dendrogram
scl_dist
scl_dendrogram
imputed_trees # imputed phylogenies to account for uncertainty
expected_assemblages # community matrix



# NULL COMMUNITIES ####

# we are going to define one global pool: any species can occur in any grid (Swenson et al. 2012)
# null.model=c('richness','independentswap'). Makes sense to compare to a similar richness structure

# number of null communities computed
Xnull=499
Nswaps=10000

# # expected assemblages
# expected_assemblages[1:5, 1:5]
# # randomize matrix and save
# null_exp_ass <- list()
# for (i in 1:Xnull) {
#   # randomize assemblage
#   null_exp_ass[[i]] <- expected_assemblages[,-c(1:3)] %>%
#     # as it goes really fast, increase Nswaps to 50000
#     picante::randomizeMatrix(null.model='richness') %>% # iterations=Nswaps
#     round(2) # needed to save space
#   # progress
#   print(paste(round(i/Xnull, 4)*100, '%'))
# }
# save(null_exp_ass, file="results/null_ass/null_exp_ass_rich.RData")



# CALCULATE ALPHA DIVERSITY ####


# calculate Faith and Pmpd as the mean values across N imputed trees to account for imputation bias
Ntrees <- 10

# site x species presence matrix to subset communities from
expected_assemblages

# relative abundances
t <- which(colnames(expected_assemblages)%in%c('cell','x','y'))
ra_expected_assemblages <- expected_assemblages[,-t]/rowSums(expected_assemblages[,-t])
table(rowSums(ra_expected_assemblages))


# results df
exp_indices <- expected_assemblages[,1:3]

# richness using occurrence probability
exp_indices$richness <- rowSums(expected_assemblages[,-c(1:3)])

# shannon diversity
exp_indices$shannon <- vegan::diversity(ra_expected_assemblages, index="shannon")

# functional composition
temp <- BAT::cwm(ra_expected_assemblages,
                 scl_traits[,c('height','blade_area','nutlet_volume')],
                 abund=TRUE, na.rm=FALSE)
exp_indices$cwm_height <- temp[,'height']
exp_indices$cwm_blade <- temp[,'blade_area']
exp_indices$cwm_nutlet <- temp[,'nutlet_volume']

# checkpoint
# write.table(exp_indices, 'results/exp_alpha_div.txt')


# life form: proportion of annual species
exp_indices$n_ann <- NA
exp_indices$n_per <- NA
v_annuals <- scl_traits %>% subset(life_form_simp==1) %>% rownames()
v_perennials <- scl_traits %>% subset(life_form_simp==2) %>% rownames()
for (i in 1:nrow(exp_indices)) {
  # retain present taxa
  temp <- ra_expected_assemblages[i,]
  # find species present
  temp <- names(temp[,colSums(temp)>0])
  # get number of annual species
  exp_indices$n_ann[i] <- length(which(temp %in% v_annuals))
  # get number of perennial species
  exp_indices$n_per[i] <- length(which(temp %in% v_perennials))
  # progress
  print(i)
}

# functional indexes
exp_indices$Frich <- BAT::alpha(ra_expected_assemblages, scl_dendrogram)[,'Richness']
exp_indices$Fmpd <- picante::mpd(ra_expected_assemblages,
                                 as.matrix(cophenetic(scl_dendrogram)),
                                 abundance.weighted=TRUE)


# checkpoint
# write.table(exp_indices, 'results/exp_alpha_div.txt')


# # phylogenetic diversity
# mat_Prich <- matrix(nrow=nrow(ra_expected_assemblages), ncol=Ntrees)
# mat_Pmpd <- mat_Prich
# for (i in 1:Ntrees) {
#   mat_Prich[,i] <- BAT::alpha(ra_expected_assemblages[,-c(1:3)], sample(imputed_trees,1)[[1]])
#   mat_Pmpd[,i] <- picante::mpd(ra_expected_assemblages[,-c(1:3)],
#                                as.matrix(cophenetic(sample(imputed_trees,1)[[1]])),
#                                abundance.weighted=TRUE)
#   print(i)
# }
# 
# # save matrices
# exp_mat_Pdiv <- list()
# exp_mat_Pdiv[['Prich']] <- mat_Prich
# exp_mat_Pdiv[['Pmpd']] <- mat_Pmpd
# save(exp_mat_Pdiv, file="results/exp_mat_Pdiv.RData")
load("results/exp_mat_Pdiv.RData")

# calculate mean and save
exp_indices$Prich <- rowMeans(mat_Prich, na.rm=F)
exp_indices$Pmpd <- rowMeans(mat_Pmpd, na.rm=F)


# checkpoint
# write.table(exp_indices, 'results/exp_alpha_div.txt')
exp_indices <- read.csv("results/exp_alpha_div.txt", sep="")

# pairs(exp_indices[,c("richness","Frich","Fmpd","Prich","Pmpd")],
#       upper.panel=NULL, main='pairs(expected alpha diversity)')



# NULL MODELS: compute indices across null communities to test for significance ####

alpha_exp_null <- matrix(nrow=nrow(expected_assemblages), ncol=Xnull) # results for every combination
rownames(alpha_exp_null) <- rownames(expected_assemblages)

null.Frich <- alpha_exp_null
null.Fmpd <- alpha_exp_null
null.Prich <- alpha_exp_null
null.Pmpd <- alpha_exp_null

# random communities via richness
load("results/null_ass/null_exp_ass_rich.RData")

# check match
length(scl_dendrogram$labels) == length(colnames(null_exp_ass[[1]]))
table(scl_dendrogram$labels %in% colnames(null_exp_ass[[1]]))
table(rownames(exp_indices)==rownames(null_exp_ass[[1]]))


for (i in 1:Xnull) {

  # relative abundances
  ra_exp_ass1 <- null_exp_ass[[i]]/rowSums(null_exp_ass[[i]])

  # Frich
  null.Frich[,i] <- BAT::alpha(ra_exp_ass1, scl_dendrogram)

  # Fmpd
  null.Fmpd[,i] <- picante::mpd(ra_exp_ass1,
                                as.matrix(cophenetic(scl_dendrogram)), abundance.weighted=TRUE)

  # temporal matrices to store phylo diversity
  mat_Prich <- matrix(nrow=nrow(ra_exp_ass1), ncol=Ntrees)
  mat_Pmpd <- matrix(nrow=nrow(ra_exp_ass1), ncol=Ntrees)

  for (t in 1:Ntrees) {
    mat_Prich[,t] <- BAT::alpha(ra_exp_ass1, sample(imputed_trees,1)[[1]])
    mat_Pmpd[,t] <- picante::mpd(ra_exp_ass1,
                                 as.matrix(cophenetic(sample(imputed_trees,1)[[1]])),
                                 abundance.weighted=TRUE)
    print(paste(t))
  }

  null.Prich[,i] <- rowMeans(mat_Prich, na.rm=TRUE)
  null.Pmpd[,i] <- rowMeans(mat_Pmpd, na.rm=TRUE)


  # progress
  print(paste('--- ', round(i/Xnull*100, 2), '% ---', sep=''))

}


# # store in list
# alpha_exp_null_list <- list() # save
# alpha_exp_null_list[['Frich']] <- null.Frich
# alpha_exp_null_list[['Fmpd']] <- null.Fmpd
# alpha_exp_null_list[['Prich']] <- null.Prich
# alpha_exp_null_list[['Pmpd']] <- null.Pmpd
# 
# # save
# save(alpha_exp_null_list, file='results/alpha_exp_null_list.Rdata')
load("results/alpha_exp_null_list.RData")



# identify significant richness and diversity patterns
exp_indices$perc_Frich <- NA
exp_indices$perc_Fmpd <- NA
exp_indices$perc_Prich <- NA
exp_indices$perc_Pmpd <- NA

for (i in 1:nrow(exp_indices)) {
  if (exp_indices$richness[i]>1) {
    exp_indices$perc_Prich[i] <- ecdf(alpha_exp_null_list$Prich[i,])(exp_indices$Prich[i])
    exp_indices$perc_Frich[i] <- ecdf(alpha_exp_null_list$Frich[i,])(exp_indices$Frich[i])
    exp_indices$perc_Fmpd[i] <- ecdf(alpha_exp_null_list$Fmpd[i,])(exp_indices$Fmpd[i])
    exp_indices$perc_Pmpd[i] <- ecdf(alpha_exp_null_list$Pmpd[i,])(exp_indices$Pmpd[i])
  }
  # progress
  print(i)
}

table(exp_indices$perc_Frich<0.05)

# # save
# write.table(exp_indices, 'results/exp_alpha_div.txt')


