


# COMPILE AND SUMMARIZE MODEL PERFORMANCE
source('SDM/scripts/SDM organize analyses.R')



# EXTRACT AUC & PLOT ####


# IDW
sum_IDW <- read.csv("SDM/results/IDW_AUC.txt", sep="")


# ENSEMBLE
sum_ENS <- list.files('SDM/results/ensemble/evaluations/ensemble', full.names=TRUE)
t <- subset(read.table(sum_ENS[[1]]), metric.eval%in%c('ROC','BOYCE'))
for (i in 2:length(sum_ENS)) { t <- rbind(t, subset(read.table(sum_ENS[[i]]), metric.eval%in%c('ROC','BOYCE'))) }
sum_ENS <- t


# MAXENT
load("SDM/results/maxent/allmodels.Rdata")
allmodels = allmodels[-which(sapply(allmodels, is.null))]
load("SDM/results/maxent/bestmodels.Rdata")
sum_MAX <- do.call(rbind, bestmodels)


# MERGE
models_AUC <- data.frame(method=c(rep('MAXENT', length(sum_MAX$auc.val.avg)),
                                  rep('IDW', length(sum_IDW$AUC)),
                                  rep('ENSEMBLE', length(sum_ENS$calibration[sum_ENS$metric.eval=='ROC']))),
                         AUC=c(sum_MAX$auc.val.avg, sum_IDW$AUC, sum_ENS$calibration[sum_ENS$metric.eval=='ROC']))

models_AUC$method <- factor(models_AUC$method, levels=c('IDW','MAXENT','ENSEMBLE'))

ggplot(aes(x=method, y=AUC), data=models_AUC) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=11, color='black'), axis.text.x=element_text(size=11, color='black'))

models_AUC %>% group_by(method) %>% summarise(mean=mean(AUC), sd=sd(AUC))



# VARIABLE IMPORTANCE ####

# I think it is not possible to compare between algorithms (MAX vs ENS) because variable importance is estimated differently

# maxent
load("SDM/results/maxent/varimp_bestmodels.Rdata")
vi_MAX <- do.call(rbind, varimp_bestmodels)
head(vi_MAX)

ggplot(aes(x=variable, y=permutation.importance), data=vi_MAX[vi_MAX$permutation.importance>0,]) +
  geom_boxplot() +
  theme_bw() +
  ylab('Variable importance MAXENT') +
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=11, color='black'),
        axis.text.x=element_text(size=11, color='black', angle=-90, vjust=0, hjust=0))


# ensemble: import single models, filter x boyce, retain variable importance, filter models
sing_mod <- list.files('SDM/results/ensemble/evaluations/simple', full.names=TRUE)
t <- read.table(sing_mod[[1]]) %>% subset(metric.eval=='BOYCE')
for (i in 2:length(sing_mod)) { t <- rbind(t, read.table(sing_mod[[i]]) %>% subset(metric.eval=='BOYCE')) }
sing_mod <- t %>% subset(validation>0.7) %>% select(full.name) %>% deframe()

vi_ENS <- list.files('SDM/results/ensemble/var_importance', full.names=TRUE)
t <- read.table(vi_ENS[[1]])
for (i in 2:length(vi_ENS)) { t <- rbind(t, read.table(vi_ENS[[i]])) }
vi_ENS <- t %>% subset(full.name %in% sing_mod)

head(vi_ENS)
ggplot(aes(x=expl.var, y=var.imp), data=vi_ENS) +
  geom_boxplot() +
  theme_bw() +
  ylab('Variable importance ENSEMBLE') +
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=11, color='black'),
        axis.text.x=element_text(size=11, color='black', angle=-90, vjust=0, hjust=0))


