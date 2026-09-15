# AUTHOR : Pascale Lemieux
# Comparison of Specificity screen with PPI screen

library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

specificity <- 
  read_csv('med_score_filtered.csv')

ppi <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ppi_signif_score.csv')

ppi$PBD <- ppi$family


specificity[specificity$pool == 'P4', 'PBD'] <- '363'
specificity[specificity$pool == 'P5', 'PBD'] <- '366'
specificity[specificity$pool == 'P6', 'PBD'] <- '385'

specificity[specificity$family == 'reference', 'PBD']<-'reference'

comp_screens <- 
merge(
  specificity[specificity$condition == 'MTX',],
  ppi,
  by.x = c('aa_seq', 'PBD', 'family'),
  by.y =  c('aa_seq', 'PBD', 'family'),
  all = F
)

specificity%>%
filter(condition == 'MTX')%>%
  inner_join(ppi, by =c('family', 'PBD', 'aa_seq'))->comp_screens

specificity%>%
  select(family, aa_seq, PBD)%>%
  filter(family=='reference')

ppi%>%
  select(family, aa_seq, PBD)%>%
  filter(family=='reference')


val_btw_screen<-
comp_screens%>%
  filter(pool %in% c('P4', 'P5', 'P6'))%>%
ggplot()+
  #facet_grid(cols = vars(family))+
  geom_abline(color = 'grey45', linetype = 'dashed')+
  geom_point(aes(med_norm, med_norm_sel_coef), color = 'darkgrey', alpha= 0.2)+
  geom_point(data =comp_screens[comp_screens$family == 'reference' & comp_screens$pool != 'P1', ], 
             aes(med_norm, med_norm_sel_coef), color = '#ff0000ff', alpha =0.6)+
  labs(x='PPI score\n(specificity assay)', y='PPI score\n(PPI assay)')+
  stat_cor(aes(med_norm, med_norm_sel_coef), method = 'spearman', cor.coef.name = 'rho', size = 3, 
           label.sep = '\n')+
  lims(x = c(-0.5, 1.2), y = c(-0.5,1.2))+
  t

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/validation_btw_screen.png', width = 4, height = 4)
saveRDS(val_btw_screen, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/validation_btw_screen.rds')

