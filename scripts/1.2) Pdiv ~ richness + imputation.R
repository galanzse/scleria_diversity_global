

# Impute Scleria phylogeny might have a great impact on phylo diversity calculation (we are imputing 107 species out of 239) or not (because we are imputing at the section levels which reduces uncertainty greatly)

# If not, we can save a lot of time by using just one phylogeny


source('scripts/1) import and prepare data.R')

library(BAT) #Frich
# devtools::install_github("iramosgutierrez/randtip")
library(randtip) # tree imputation
library(ggpubr)



# Lets impute 100 trees ####

# Ntrees = 100
# 
# scl_phylogeny <- data_final[['phylogeny']]
# 
# # used simplified epithets
# scl_phylogeny$tip.label <- scl_taxonomy$epithet[match(scl_phylogeny$tip.label, scl_taxonomy$scientific_name)]
# 
# scl_phylogeny <- drop.tip(scl_phylogeny, scl_phylogeny$tip.label[which(!(scl_phylogeny$tip.label %in% scl_occurrences$epithet))]) # drop species without occurrences
# 
# scl_taxonomy <- scl_taxonomy %>% filter(epithet %in% v_species) # infrageneric classification, treat sections as genera for imputation
# scl_taxonomy$sectxspp <- paste(scl_taxonomy$section, scl_taxonomy$epithet)
# 
# phylogeny_sect <- scl_phylogeny # change names in tree
# intree <- scl_taxonomy[scl_taxonomy$epithet %in% phylogeny_sect$tip.label,]
# phylogeny_sect$tip.label <- intree$sectxspp[order(match(intree$epithet, phylogeny_sect$tip.label))]
# 
# back.tree <- phylogeny_sect # backbone phylogeny
# class(back.tree)
# is.ultrametric(back.tree)
# 
# sp.list <- scl_taxonomy$sectxspp # species list
# 
# # build info df
# my.info.noranks <- randtip::build_info(species=sp.list, tree=back.tree, mode='backbone', find.ranks=FALSE)
# 
# # check: yields errors but it is ok
# my.check <- randtip::check_info(my.info.noranks, tree=back.tree)
# my.check$Typo[my.check$taxon=='Hymenolytrum_ramosa'] <- FALSE
# my.check$Typo.names[my.check$taxon=='Hymenolytrum_ramosa'] <- NA
# 
# # define ranks
# my.input.noranks <- randtip::info2input(my.info.noranks, back.tree)
# 
# # to change to scientific names after imputation
# scl_taxonomy$sectxspp <- gsub(' ', '_', scl_taxonomy$sectxspp)
# 
# imputed_trees <- list()
# for (i in 1:Ntrees) {
#   imp_tree <- randtip::rand_tip(input=my.input.noranks, tree=back.tree,
#                                 rand.type='random', # default
#                                 respect.mono=T, # all our clades are monophyletic so not relevant
#                                 prob=F, # branch selection probability is equiprobable
#                                 use.stem=T, # the stem branch can be considered as candidate for binding
#                                 prune=F,
#                                 verbose=F)
# 
#   # (!!!) Nnode(imp_tree) == length(imp_tree$tip.label)
# 
#   imp_tree$tip.label <- scl_taxonomy$epithet[order(match(scl_taxonomy$sectxspp, imp_tree$tip.label))] # retrieve original names
# 
#   imputed_trees[[i]] <- imp_tree # save in list
#   ck_n <- deframe(table(imp_tree$tip.label%in%scl_taxonomy$epithet))
#   print(paste('x',i,'   ', 'Ntip=', ck_n, sep='')) # progress
# }
# 
# 
# rm(imp_tree, ck_n, my.info.noranks, my.check, my.input.noranks, phylogeny_sect, intree)
# 
# save(imputed_trees, file="results/imputed_trees.RData")
load("results/imputed_trees.RData")



# Relationship of mean and mad indexes with richness ####

# Lets impute N communities of R richness levels I times 
R=20
N=R*10
I=99

# randomly fill each row with the R level we want 
comm_pool <- observed_assemblages[1:N,-c(1:3)]
comm_pool[] <- 0
rl <- rep(1:R, 10) # length(rl)==nrow(comm_pool)

for (i in 1:nrow(comm_pool)) {
  # randomly select positions according to richness
  pos <- sample(1:ncol(comm_pool), rl[i])    
  # fill those positions with 1 (presence)
  comm_pool[i,pos] <- 1
}
# table(rowSums(comm_pool)==rl)


# calculate phylogenetic richness and diversity for each community I times
imp_importance <- data.frame(richness=rowSums(comm_pool),
                             Prich=NA, mad_Prich=NA,
                             Pmpd=NA, mad_Pmpd=NA)

for (i in 1:N) {
  
  # vector to store indexes after each imputation
  v_Prich <- vector(length=I)
  v_Pmpd <- vector(length=I)
  
  # extract I trees from list of imputed phylogenies
  rdm_trees <- sample(imputed_trees, I, replace=F)
  
  # calculate indexes
  for (i2 in 1:I) {
    v_Prich[i2] <- BAT::alpha(comm_pool[i,], rdm_trees[[i2]])[,'Richness']
    v_Pmpd[i2] <- picante::mpd(comm_pool[i,], as.matrix(cophenetic(rdm_trees[[i2]])), abundance.weighted=F)
  }
  
  # compute median and median absolute deviation
  imp_importance$Prich[i] <- median(v_Prich)
  imp_importance$mad_Prich[i] <- mad(v_Prich)
  imp_importance$Pmpd[i] <- median(v_Pmpd)
  imp_importance$mad_Pmpd[i] <- mad(v_Pmpd)
  
  # progress
  print(paste(round(i/N*100, 2), '%'))

}

# save
# write.table(imp_importance, 'results/imputation_importance.txt')


# Lets explore the data:
a1<- ggplot(aes(x=richness, y=Prich), data=imp_importance) +
  geom_point() +
  geom_errorbar(aes(ymin=Prich-mad_Prich, ymax=Prich+mad_Prich)) +
  theme_bw()
a2 <- ggplot(aes(x=richness, y=Pmpd), data=imp_importance) +
  geom_point() +
  geom_errorbar(aes(ymin=Pmpd-mad_Pmpd, ymax=Pmpd+mad_Pmpd)) +
  theme_bw()
a3<- ggplot(aes(x=Prich, y=mad_Prich), data=imp_importance) +
  geom_point() +
  xlab('Prich') + ylab('median absolute deviation') +
  theme_bw()
a4 <- ggplot(aes(x=Pmpd, y=mad_Pmpd), data=imp_importance) +
  geom_point() +
  xlab('Pmpd') + ylab('median absolute deviation') +
  theme_bw()

ggarrange(a1, a2, a3, a4, ncol=2, nrow=2)


# Imputation is not that important


