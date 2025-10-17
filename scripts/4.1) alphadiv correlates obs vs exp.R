

library(tidyverse)


# data
obs_alphadiv <- read.csv("results/test_obs_alphadiv_bioclim.csv", row.names=1)
obs_alphadiv$type <- 'observed'
obs_alphadiv$P_value[obs_alphadiv$P_value<0.05] <- 'sign.'
obs_alphadiv$P_value[obs_alphadiv$P_value!=('sign.')] <- 'n.s.'

exp_alphadiv <- read.csv("results/test_exp_alphadiv_bioclim.csv", row.names=1)
exp_alphadiv$type <- 'expected'
exp_alphadiv$P_value[exp_alphadiv$P_value<0.05] <- 'sign.'
exp_alphadiv$P_value[exp_alphadiv$P_value!=('sign.')] <- 'n.s.'


# compare correlation coefficients and pvals
corr_vals <- rbind(obs_alphadiv, exp_alphadiv) %>%
  dplyr::select(ALPHA, BIOCLIM, corr, type) %>%
  pivot_wider(names_from=type, values_from=corr) %>%
  na.omit()

# do significancy match?
corr_vals$P_sim <- NA
for (i in 1:nrow(corr_vals)) {
  t1 <- exp_alphadiv$P_value[exp_alphadiv$ALPHA==corr_vals$ALPHA[i] & exp_alphadiv$BIOCLIM==corr_vals$BIOCLIM[i]]
  t2 <- obs_alphadiv$P_value[obs_alphadiv$ALPHA==corr_vals$ALPHA[i] & obs_alphadiv$BIOCLIM==corr_vals$BIOCLIM[i]]
  if (t1==t2) { corr_vals$P_sim[i] <- 'Match' } else { corr_vals$P_sim[i] <- 'Do not match' }
}
table(corr_vals$P_sim)/sum(table(corr_vals$P_sim))

corr_vals$ALPHA <- as.factor(corr_vals$ALPHA)
levels(corr_vals$ALPHA) <- c("CWM Blade area", "CWM Height", "CWM Nutlet volume", "Fmpd", "Frich", "Pmpd", "Prich", "Richness")

ggplot(aes(x=observed, y=expected, color=ALPHA), data=corr_vals) +
  geom_point(aes(shape=P_sim)) +
  scale_shape_manual(values=c(4,16)) +
  theme_classic() +
  xlim(-0.6, 0.6) + ylim(-0.8, 0.8) +
  xlab('observed correlation') + ylab('expected correlation') +
  # geom_smooth(method='lm', se=F) +
  geom_abline(intercept=0, slope=1, linetype=1, size=1.5) +
  geom_hline(yintercept=0, linetype=2) + geom_vline(xintercept=0, linetype=2) +
  labs(shape="P value", color="Index")


