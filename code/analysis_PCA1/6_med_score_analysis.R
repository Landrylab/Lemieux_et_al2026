# AUTHOR : Pascale Lemieux
### Compare PPI and Avail score to stickiness and PWM score ###
#import packages

library(tidyverse)
library(Biostrings)
library(magrittr)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(jsonlite)
library("ggsci")

source('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/code/functions.R')

# set working directory

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/data/PCA1/')

# import selection coefficient
med_ppi <- read_csv('med_ppi_score.csv')
med_avail <- read_csv('med_avail_score.csv')
ppi_signif <- read_csv('ppi_signif_score.csv')

# merge by aa sequence

comp_ppi_avail <- 
merge(med_ppi[, -c(2:8)], 
      med_avail, 
      by = 'aa_seq', 
      suffixes = c('.ppi', '.avail'))

comp_ppi_avail$stop <- 
  c(comp_ppi_avail$deg1 == '*' | comp_ppi_avail$deg2 == '*')

comp_ppi_avail$max_pep <- 
(comp_ppi_avail$deg1 == comp_ppi_avail$max1 & comp_ppi_avail$deg2 == comp_ppi_avail$max2)


comp_ppi_avail$name <- 
  gsub(comp_ppi_avail$family, pattern = 'PRM_0', replacement = '')


#Fig3D<- 
ggplot(comp_ppi_avail)+
  geom_point(data = comp_ppi_avail, 
            aes(med_norm.avail, med_norm.ppi, color = name, shape = stop, size = stop, alpha=stop))+
  scale_color_viridis_d(option='G')+
  scale_size_manual(values = c(1, 2))+
  scale_alpha_manual(values = c(0.4, 1))+
  scale_shape_manual(values = c(1, 4))+
  ylab('PPI score')+
  xlab('Avail. score')+
  guides(color = guide_legend(title = 'peptide family', position = 'bottom'), 
         shape = guide_legend(nrow = 2, title ='stop codon', position = 'bottom'), 
         size = 'none', 
         alpha = 'none')+
 t


#saveRDS(Fig3D, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/FigAvailD.rds')


# Compare availability with stickyness

sticky <- 
  read_table('../sticky_score.txt')

sub_comp <- 
  med_avail[!is.na(med_avail$med_norm), ]

aa_coding <- 
  strsplit(unlist(sub_comp$aa_seq), split = '*', fixed = T)

aa_seq <- lapply(aa_coding, `[[`, 1)

aa <- lapply(aa_seq, strsplit, split = '')

# stop peptide with first stop codon
#sub_coding <- comp_ppi_avail[!(comp_ppi_avail$stop | is.na(comp_ppi_avail$stop)), ]

for (i in 1:length(aa)) {
  med_avail[i, 'sticky_score'] <- get_stickyness(aa[[i]], sticky)
}

med_avail$name <- 
  gsub(med_avail$family, pattern = 'PRM_0', replacement = '')

med_avail$stop <- 
  grepl(med_avail$aa_seq, pattern ='*', fixed = T)

# absolute stickyness
#Fig3E <- 
ggplot(med_avail)+
  #facet_wrap(vars(family))+
  geom_point(data = med_avail[!med_avail$stop, ],
             aes(x = med_norm, sticky_score, color = name), shape = 1, alpha = 0.5, size = 1)+
  geom_point(data = med_avail[med_avail$stop, ],
             aes(x = med_norm, sticky_score, color = name), shape = 4, alpha = 1, size = 2)+
  labs(x = 'Avail. score', y = 'absolute stickyness')+
  scale_color_viridis_d(option = 'G')+
  #scale_color_frontiers()+
  scale_shape_manual(values= c(20,4))+
  stat_cor(aes(x = med_norm, sticky_score), 
           method = 'spearman', label.sep = '\n', cor.coef.name = 'rho', label.x = 0.8, size = 3)+
  t+
  #theme(legend.position = 'none')+
  guides(color= 'none', 
         shape = 'none')


saveRDS(Fig3E, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigAvail_vsSticky.rds')

# compute delta stickyness
med_avail%>%
  dplyr::filter(deg1 == max1& deg2 == max2 & family != 'reference')-> reference_avail

colnames(reference_avail)[c(9, 13)] <- c('med_norm.ref', 'sticky_score.ref')

med_avail <- 
left_join(med_avail, 
          reference_avail[,c(8,9,13)], 
          by = 'family')

med_avail%>%
  mutate(delta_avail = med_norm - med_norm.ref, delta_sticky = sticky_score - sticky_score.ref)->delta_avail


# delta stickyness vs delta avail
ggplot(delta_avail)+
  #facet_wrap(vars(family))+
  geom_point(aes(x = med_norm, delta_sticky, color = family), alpha = 0.5)+
  stat_cor(aes(x = med_norm, delta_sticky))+
  theme_classic2()


# vizualise stickiness between max peptide and peptide with stronger ppi
ppi_signif_stronger <- 
ppi_signif[ppi_signif$binding == 'stronger', ]

comp_max_stronger <- 
merge(med_avail, 
      ppi_signif_stronger[, -c(2:7)], 
      by = 'aa_seq', 
      all.y =T, 
      all.x =F)


comp_max_stronger%>%
  dplyr::group_by(name)%>%
  dplyr::mutate(med_sticky = median(sticky_score, na.rm =T))%>%
  dplyr::select(name, sticky_score.ref, med_sticky)%>%
  na.omit()%>%
  unique()-> med_sticky


Fig3F <- 
ggplot(med_sticky) +
  geom_segment(aes(x = sticky_score.ref, xend = med_sticky,
                   y = name, yend = name)) +
  geom_point(aes(x = sticky_score.ref, y = name), color = 'black', size = 3.6)+ 
  geom_point(aes(x = sticky_score.ref, y = name, color = name), size = 3) +
  scale_color_manual(values = c( "#2E1E3CFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  geom_point(aes(x = med_sticky, y = name), size = 3, shape = 4)+
  ylab('peptide family')+
  xlab('absolute stickiness')+
  t+guides(color = 'none')

saveRDS(Fig3F, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig_maxPPI_vs_sticky.rds')


## comparison with PWM

# Import PMW scores
PWM <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/data/PWM_updated.rds')

signif_data <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ppi_signif_score.csv')
med_norm<-read_csv('med_ppi_score.csv')

assay_PWM<-PWM[names(PWM) %in% unique(med_norm$family)]

# create object to save the PWM scores
all_score <- tibble()

# loop through the PBD ID and compute PWM score of all peptide against each PWM
fam<-unique(med_norm$family)[-5]

for(n in fam){
  
  sub_PWM <- 
    assay_PWM[names(assay_PWM) == n][[1]]
  
  med_norm%>%
    filter(family == n)%>%
    select(aa_seq)%>%
    unlist()%>%
    unique()->seq_pep
  
  score<-vector('numeric', length(seq_pep))
  
  for(j in 1:length(seq_pep)){
    get_PWM_from_seq(seq_pep[j], sub_PWM = sub_PWM) -> score[j]
  }
  
  x<-bind_cols(seq_pep, unlist(score), n)
  
  all_score <- bind_rows(all_score, x)
}

colnames(all_score)<- c('aa_seq', 'PWM_score', 'PBD')



all_score%>%
  right_join(med_norm, join_by(aa_seq, PBD==family))->med_norm


signif_data%>%
  select(aa_seq, binding)%>%
  right_join(med_norm)->med_norm


med_norm%>%
  filter(max)%>%
  select(PBD, PWM_score)->ref_PWM

colnames(ref_PWM)[2]<-'PWM_score_max'

med_norm%<>%
  left_join(ref_PWM, join_by(PBD))%>%
  mutate(rel_PWM_score = PWM_score/PWM_score_max)

med_norm$PBD<-
  gsub('PRM_0', 'PBD ', med_norm$PBD)


med_norm%>%
  filter(PBD !='reference')%>%
  filter(PBD !='PBD 366')%>%
  ggplot()+
  facet_wrap(vars(PBD))+
  #geom_smooth(aes(x= rel_PWM_score, y = med_norm), color='black', fill = 'grey75', method = 'loess')+
  geom_hline(data = med_norm[med_norm$max & !(med_norm$PBD %in% c('reference', 'PBD 366')), ],
             aes(yintercept = med_norm), 
             color = 'grey45', linetype = 'dashed')+
  geom_point(aes(x= rel_PWM_score, y = med_norm, color = binding, alpha = binding, shape = binding))+
  scale_shape_manual(values = c(1, 19, 1))+
  scale_color_manual(values = c('grey','#4e67c8ff' , '#373737' ))+
  scale_alpha_manual(values =c(0.7, 1, 1))+
  stat_cor(aes(x= PWM_score, y = med_norm), method = 'spearman', size = 4, label.sep = '\n', 
           cor.coef.name = 'rho', label.x = 0.75 , label.y =0.8)+
  labs(x = 'PWM score', y = 'PPI score')+
  # lims(x = c(70,140), y = c(0, 1))+
  t+
  #theme(legend.position = 'none')+
  guides(color = guide_legend(nrow = 1, position = 'bottom', title = ''), 
         shape = 'none', 
         alpha = 'none')-> figS_comp_PWM

write_csv(med_norm, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/PWM_scores.csv')
med_norm<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/PWM_scores.csv')

saveRDS(figS_comp_PWM, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Supp_comp_with_PWM_exp1.rds')

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS5_PWM_exp1.svg', figS_comp_PWM,
       width = 8, height = 5)

med_norm%>%
  filter(PBD =='PBD 366')%>%
  ggplot()+
  facet_wrap(vars(PBD))+
  geom_smooth(aes(x= rel_PWM_score, y = med_norm), color='black', fill = 'grey75', method = 'loess')+
  geom_hline(data = med_norm[med_norm$max & (med_norm$PBD %in% c('366')), ],
             aes(yintercept = med_norm), 
             color = 'grey45', linetype = 'dashed')+
  geom_point(aes(x= rel_PWM_score, y = med_norm, color = binding, alpha = binding, shape = binding))+
  scale_shape_manual(values = c(1, 19, 1))+
  scale_color_manual(values = c('grey','#4e67c8ff' , '#373737'))+
  scale_alpha_manual(values =c(0.7, 1, 1))+
  stat_cor(aes(x= PWM_score, y = med_norm), method = 'spearman', size = 4, label.sep = '\n', 
           cor.coef.name = 'rho', label.x = 0.95 , label.y =0.3)+
  labs(x = 'PWM score', y = 'PPI score')+
  #lims(x = c(0.85,1), y = c(0, 1))+
  t+
  #theme(legend.position = 'none')+
  guides(color = guide_legend(nrow = 1, position = 'bottom', title = ''), 
         shape = 'none', 
         alpha = 'none')-> Fig2G

saveRDS(Fig2G, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2G.rds')

med_avail%>%
  dplyr::select(aa_seq, med_norm)%>%
  right_join(med_norm, join_by(aa_seq), suffix = c('.ppi', '.avail'))->test_model


avail_PWM_model<-lm(med_norm.ppi ~ rel_PWM_score + med_norm.avail, data = test_model)

PWM_model<-lm(med_norm.ppi ~ rel_PWM_score, data = test_model)

summary(avail_PWM_model)

summary(PWM_model)
