# AUTHOR : Pascale Lemieux
### Comparison of PCA scores with position weight matrices obtained by phage display ###

#import packages and custom functions

library(tidyverse)
library(Biostrings)
library(magrittr)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(jsonlite)
library("ggsci")
library(magrittr)
library(GGally)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# import median scores
signif <- read_csv('signif_score_welch20_filter.csv')


# Import PMW scores
PWM <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/data/PWM_updated.rds')


# keep only peptide families of interest from the PWM dataset
fam<-str_c('PRM_0', unlist(unique(signif[signif$family !='reference', 'family'])))
assay_PWM <- 
  PWM[names(PWM) %in% fam]

# format PCA score dataset
pep_deg <- 
  signif[, c('aa_seq', 'pos1', 'pos2', 'pos','family', 'signif')]

pep_deg%<>%
  mutate(family = str_c('PRM_0', family))

pep_deg %<>%mutate_at(c(2:4), as.integer)

# create object to sae the PWM scores
all_score <- tibble()

# compute PWM scores only for the peptide families included in the specificity assay
seq_pdz<-
unlist(unique(signif[signif$pool %in% c('P4', 'P5', 'P6'), 'aa_seq' ]))

# loop through the PBD ID and compute PWM score of all peptide against each PWM
for(n in c('PRM_0363', 'PRM_0366', 'PRM_0385')){
  
  sub_PWM <- 
    assay_PWM[names(assay_PWM) == n][[1]]
  
  score<-vector('numeric', length(seq_pdz))
  
  for(j in 1:length(seq_pdz)){
    get_PWM_from_seq(seq_pdz[j], sub_PWM = sub_PWM) -> score[j]
  }
  
  x<-bind_cols(seq_pdz, unlist(score), n)
 
  all_score <- bind_rows(all_score, x)
}

colnames(all_score)<- c('aa_seq', 'score', 'PBD')

#associate PCA scores with PWM scores
all_score%<>%
mutate(
        pool = case_when(
            PBD == "PRM_0363" ~  'P4',
            PBD == "PRM_0366" ~ 'P5',
            PBD == "PRM_0385" ~ 'P6'
        )
    ) 

comp_ppi_avail <- 
left_join(all_score[, -c(3)], 
          signif[, c('aa_seq', 'pool','med_norm', 'peptide_max', 'assay', 'condition', 'side', 'signif', 'family')], 
          by = c('aa_seq', 'pool'))


# select only PPI scores
comp_ppi_avail%>%
filter(assay == 'PPI' & condition == 'MTX') %>%
filter(pool %in% c('P4', 'P5', 'P6') & family %in% c(363, 366, 385))%>%
unique()-> comp_PWM


comp_PWM%<>%
  group_by(pool)%>%
  mutate(PWM_score_max = max(score, na.rm = TRUE))%>%
  ungroup()%>%
  mutate(rel_PWM_score = score/PWM_score_max)

write_csv(comp_PWM, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/PWM_score.csv')


# plot PPI scores vs PWM scores 
comp_PWM%>%
  filter(family == 363 & pool == 'P4' | family == 366 & pool == 'P5' | family == 385 & pool == 'P6')->sub_data
  
comp_PWM%>%
ggplot(aes(x= rel_PWM_score, y = med_norm))+
  stat_density_2d(aes(fill =..level..),
    contour_var = "ndensity",
    geom = "polygon", color = 'white', linewidth = 0.2)+
  scale_fill_continuous(palette = c('white', 'grey30', 'black'))+
  stat_cor(aes(x= rel_PWM_score, y = med_norm), method = 'spearman', size = 3, label.sep = '\n', 
  cor.coef.name = 'rho', label.x = 0.35 , label.y =0.85)+
  labs(x = 'PWM score', y = 'PPI score')+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.1))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.15, 1.1))+
  t+
  #theme(legend.position = 'none')+
  guides(fill = guide_colorbar(position = 'bottom', title = 'density')) -> PPI_all_vs_PWM

saveRDS(PPI_all_vs_PWM, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ppi_all_vs_PWM.rds')

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ppi_all_vs_PWM.png',
 width = 3.5, height = 3.6)

#correlation plot between the PBD-peptide family combination

#pivot wider
comp_PWM$PBD<-comp_PWM$pool

comp_PWM$PBD<-gsub(pattern = 'P4', replacement = 363, comp_PWM$PBD)
comp_PWM$PBD<-gsub(pattern = 'P5', replacement = 366, comp_PWM$PBD)
comp_PWM$PBD<-gsub(pattern = 'P6', replacement = 385, comp_PWM$PBD)

comp_PWM%>%
  select(PBD, family, aa_seq, score, med_norm)%>%
  pivot_wider(names_from = PBD, values_from = c(score, med_norm))%>%
  na.omit()->wide_comp

comp_PWM$family<-
as.factor(comp_PWM$family)

comp_PWM$PBD<-
  as.factor(comp_PWM$PBD)

comp_PWM%>%
  #filter(pool == 'P5')%>%
  ggcorr(columns =  c(2,4,10,11), ggplot2::aes(color = family, alpha = PBD), method = 'spearman')


# select the highest scoring peptides for each library associated to its PBD

# select the top 5 significant peptides per PBD ordered by their PPI score
comp_PWM %>%
filter(family == '363' & pool == 'P4'|
        family == '366' & pool == 'P5'|
        family == '385' & pool == 'P6') %>%
  dplyr::group_by(family, pool) %>%
  filter(side == 'stronger') %>%
  slice_max(order_by = med_norm, n = 5) %>%
  ungroup() -> peptide_max_PWM

# select the top 5 peptides with the strongest PPI score (>0.9, not necessarily significant) per PBD ordered 
comp_PWM %>%
  filter(family == '363' & pool == 'P4'|
         family == '366' & pool == 'P5'|
         family == '385' & pool == 'P6') %>%
  dplyr::group_by(family, pool) %>%
  filter(med_norm>0.9) %>%
  slice_max(order_by = med_norm, n = 5) %>%
  ungroup()-> peptide_max_PWM_2

# combine peptides selected by both criteria
bind_rows(peptide_max_PWM, peptide_max_PWM_2)%>%
  unique()-> peptide_max_PWM_all

# select information for those selected peptides
comp_PWM[comp_PWM$aa_seq %in% peptide_max_PWM_all$aa_seq, ]->
  validation_PWM


# save the peptide to use for validation, the top - 5 significativly stronger binders and 
# all peptide above 0.9 PPI score
write_csv(peptide_max_PWM_all, '~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/peptide_validation_max.csv')

# Select peptide with stronger PPI but for the non design PBD

signif%>%
filter(!grepl('*', aa_seq, fixed =TRUE))%>%
filter(condition == 'MTX' & assay == 'PPI')%>%
filter(family == '363' & pool != 'P4'|
         family == '366' & pool != 'P5'|
         family == '385' & pool != 'P6')%>%
filter(med_norm> 0.46 & side =='stronger')%>%
dplyr::group_by(pool)%>%
slice_max(order_by = med_norm, n = 5) %>%
ungroup()-> pep_crossspecificity

write_csv(pep_crossspecificity, '~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/peptide_validation_specificity.csv')


### Compare pools between each other
rm(list = ls())
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

signif_data <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/signif_score_welch20_filter.csv')

val_spe<-read_csv('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/peptide_validation_specificity.csv')
val_stop<-read_csv('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/peptide_validation_stop.csv')
val_max<-read_csv('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/peptide_validation_max.csv')

val_spe$type<-'specificity'
val_stop$type<-'stop'
val_max$type<-'max'

val<-
  bind_rows(list(val_spe, val_stop, val_max), .id = 'validation')%>%
  select( aa_seq, family, pool, med_norm, type, validation)

val$type<-
  factor(val$validation, 
         levels = c('1','2','3'), 
         labels = c('specificity', 'stop', 'max'))

signif_data%>%
  filter(condition == 'MTX' & pool %in% c('P4', 'P5', 'P6'))%>%
  select(aa_seq, pool, family, med_norm, side)%>%
  dplyr::group_by(aa_seq, family)%>%
  pivot_wider(values_from = c(med_norm, side), 
              names_from = pool, values_fn = unique)-> comp_pool

signif_data%>%
  filter(condition == 'MTX' & pool %in% c('P4', 'P5', 'P6'))%>%
  filter(peptide_max)%>%
  select(aa_seq, pool, family, med_norm, side, peptide_max)%>%
  unique()->max_pep

comp_pool%<>%
  left_join(val[, c('type', 'aa_seq')], by = 'aa_seq')


comp_pool$type<-
  replace_na(as.character(comp_pool$type), 'non validation')

comp_pool$type<-
  factor(comp_pool$type, 
         levels = c('non validation', 'max', 'specificity', 'stop'))

comp_pool$validation <- comp_pool$type != 'non validation'


max_pep%<>%
  pivot_wider(names_from = pool, values_from = med_norm)


comp_pool%>%
  filter(family %in% c(363, 366, 385))%>%
  ggplot()+
  stat_density_2d(aes(x = med_norm_P4, y = med_norm_P5, color =as.factor(family)),
                  contour_var = "ndensity",
                  geom = "polygon", linewidth = 0.3, fill ='transparent')+
  geom_point(aes(x = med_norm_P4, y = med_norm_P5, color =as.factor(family), 
                 alpha = validation), size = 2, shape= 4)+
  geom_point(data = max_pep, aes(P4, P5), size =3, color = 'black')+
  geom_point(data = max_pep, aes(P4, P5, color = family), size =2.5)+
  #stat_cor(aes(x = med_norm_P4, y = med_norm_P5), method ='spearman', cor.coef.name = 'rho', 
   #        label.sep = '\n',label.x = 0.7, label.y = 0.9, size =3)+
  scale_color_manual(values = c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey'))+
  scale_alpha_manual(values = c(0, 1))+
  #scale_shape_manual(values = c(1,16,17,4))+
  labs(x = '363', y = '366')+
  lims(x = c(-0.2, 1.05), y = c(-0.2,1.05))+
  t+
  theme(legend.position ='none')->P4P5

comp_pool%>%
  filter(family %in% c(363, 366, 385))%>%
  ggplot()+
  stat_density_2d(aes(x = med_norm_P4, y = med_norm_P6, color =as.factor(family)),
                  contour_var = "ndensity",
                  geom = "polygon", linewidth = 0.3, fill ='transparent')+
  geom_point(aes(x = med_norm_P4, y = med_norm_P6, color =as.factor(family), 
                 alpha = validation), size = 2,shape= 4)+
  geom_point(data = max_pep, aes(P4, P6), size =3, color = 'black')+
  geom_point(data = max_pep, aes(P4, P6, color = family), size =2.5)+
  #stat_cor(aes(x = med_norm_P4, y = med_norm_P6), method ='spearman',  cor.coef.name = 'rho', 
   #        label.sep = '\n',label.x = 0.7, label.y = 0.9, size = 3)+
  scale_color_manual(values = c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey'))+
  scale_alpha_manual(values = c(0, 1))+
  #scale_shape_manual(values = c(1,16,17,4))+
  labs(x = '363', y = '385')+
  lims(x = c(-0.2, 1.05), y = c(-0.2,1.05))+
  t+
  theme(legend.position ='none')->P4P6


comp_pool%>%
  filter(family %in% c(363, 366, 385))%>%
  ggplot()+
  stat_density_2d(aes(x = med_norm_P5, y = med_norm_P6, color =as.factor(family)),
                  contour_var = "ndensity",
                  geom = "polygon", linewidth = 0.3, fill ='transparent')+
  geom_point(aes(x = med_norm_P5, y = med_norm_P6, color =as.factor(family), 
                 alpha = validation), size = 2, , shape= 4)+
  geom_point(data = max_pep, aes(P5, P6), size =3, color = 'black')+
  geom_point(data = max_pep, aes(P5, P6, color = family), size =2.5)+
  labs(x = '366', y = '385')+
  scale_color_manual(values = c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey'))+
  #stat_cor(aes(x = med_norm_P5, y = med_norm_P6), method ='spearman',  cor.coef.name = 'rho', 
   #        label.sep = '\n', label.x = 0.7, label.y = 0.9, size =3)+
  scale_alpha_manual(values = c(0, 1))+
  #scale_shape_manual(values = c(1,16,17,4))+
  lims(x = c(-0.2, 1.05), y = c(-0.2,1.05))+
  t+
  guides(alpha = 'none', 
         color = guide_legend(title = 'family'), 
         shape = guide_legend(title = 'peptide'))+
  theme(legend.position ='none')->P5P6


saveRDS(P4P5, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4P5.rds')
saveRDS(P4P6, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4P6.rds')
saveRDS(P5P6, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P5P6.rds')

ggsave(P4P5, filename ='~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4P5.svg', height = 4, width = 4.1)

ggsave(P4P6, filename ='~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4P6.svg', height = 4, width = 4.1)

ggsave(P5P6, filename ='~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P5P6.svg', height = 4, width = 4.1)


signif_data%<>%
  dplyr::mutate(paired = case_when((pool == 'P4' & family == '363' | 
                                     pool == 'P5' & family == '366' | 
                                     pool == 'P6' & family == '385') ~ TRUE  , 
                                    (pool == 'P4' & family != '363' | 
                                      pool == 'P5' & family != '366' | 
                                      pool == 'P6' & family != '385') ~ FALSE), 
                double_mutant = case_when(is.na(max)~TRUE,
                                          !is.na(max)~ FALSE))

# select peptide for in silico modelling
signif_data%>%
  dplyr::filter(paired & condition == 'MTX')%>%
  ggplot()+
  geom_violin(aes(med_norm, family))+
  geom_point(data = signif_data[signif_data$peptide_max & signif_data$paired &  signif_data$condition == 'MTX', ], 
             aes(med_norm, family))


signif_data%>%
  dplyr::filter(paired & condition == 'MTX')%>%
  dplyr::group_by(family)%>%
  dplyr::summarise(
    q1 = quantile(med_norm, 0.01, na.rm = TRUE),
    q10 = quantile(med_norm, 0.1, na.rm = TRUE),
    q20 = quantile(med_norm, 0.2, na.rm = TRUE),
    q30 = quantile(med_norm, 0.3, na.rm = TRUE),
    q40 = quantile(med_norm, 0.4, na.rm = TRUE),
    q50 = quantile(med_norm, 0.5, na.rm = TRUE),
    q60 = quantile(med_norm, 0.6, na.rm = TRUE),
    q70 = quantile(med_norm, 0.7, na.rm = TRUE),
    q80 = quantile(med_norm, 0.8, na.rm = TRUE),
    q90 = quantile(med_norm, 0.9, na.rm = TRUE),
    q100 = quantile(med_norm, 1, na.rm = TRUE),
    .groups = "drop"
  )-> quantile_family


peptide_select<-tibble::tibble()
q_test<-c('q1', 'q10', 'q20', 'q30', 'q40', 'q50', 'q60', 'q70', 'q80', 'q90', 'q100')

for (f in c('363', '366', '385')) {
  
  signif_data%>%
    dplyr::filter(paired & condition == 'MTX' & family ==f)%>%
    dplyr::filter(!grepl('*', aa_seq, fixed = TRUE))->sub_signif

    for (q in 1:10) {
      
      sub_q<-q_test[c(q, q+1)]
      
       quantile_family%>%
          dplyr::filter(family == f)%>%
          dplyr::select(all_of(sub_q))%>%
          unlist()->sub_quantile
      
      sub_signif%>%
        dplyr::filter(med_norm<sub_quantile[2] &  med_norm>=sub_quantile[1] & side != 'no significant difference' & double_mutant)%>%
        slice_max(order_by = med_norm, n=1, with_ties = FALSE)%>%
        select(family, aa_seq, med_norm, double_mutant, side)->q_peptide_d
      
      sub_signif%>%
        dplyr::filter(med_norm<sub_quantile[2] &  med_norm>=sub_quantile[1] & side != 'no significant difference' & !double_mutant)%>%
        slice_max(order_by = med_norm, n=1, with_ties = FALSE)%>%
        select(family, aa_seq, med_norm, double_mutant, side)->q_peptide_s
      
      if(nrow(q_peptide_d) == 0){
        
      sub_signif%>%
        dplyr::filter(med_norm<sub_quantile[2] &  med_norm>=sub_quantile[1] & double_mutant)%>%
          slice_max(order_by = med_norm, n=1, with_ties = FALSE)%>%
          select(family, aa_seq, med_norm, double_mutant, side)->q_peptide_d
      }
      
      if(nrow(q_peptide_s) == 0){
  
        sub_signif%>%
          dplyr::filter(med_norm<sub_quantile[2] &  med_norm>=sub_quantile[1]& !double_mutant)%>%
          slice_max(order_by = med_norm, n=1, with_ties = FALSE)%>%
          select(family, aa_seq, med_norm, double_mutant, side)->q_peptide_s
      }
      
        
      q_peptide_d$q<-q
      q_peptide_s$q<-q
      
      peptide_select<-bind_rows(peptide_select, q_peptide_d, q_peptide_s)
        
      
    }
  
}

signif_data%>%
  dplyr::filter(paired & condition == 'MTX')%>%
  dplyr::filter(peptide_max)%>%
  select(family, aa_seq, med_norm, side)%>%
  unique()->ref_pep

ref_pep$side<-'reference'

signif_data%>%
  dplyr::filter(paired & condition == 'MTX')%>%
  ggplot()+
  geom_violin(aes(med_norm, family))+
  geom_point(data =peptide_select, 
             aes(med_norm, family, color =side))+
  geom_point(data =ref_pep, 
             aes(med_norm, family), shape = 4)+
  scale_color_manual(values = c('red', 'grey', 'blue'), 
                     limits = c('stronger', 'no significant difference', 'weaker'))+
  theme_classic2()+
  xlab('PPI score')+
  theme(legend.position = 'bottom')



bind_rows(peptide_select, ref_pep)->xavier_peptide

colnames(xavier_peptide)[c(3,6)]<-c('PPI_score', 'q10')

write_csv(xavier_peptide, '~/PL_projects/PL_papers/PPI_optimization_paper/data/validation/peptide_docking_validation.csv')



