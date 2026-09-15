# AUTHOR : Pascale Lemieux
# visualize loss of specificity example
library(tidyverse)
library(Biostrings)
library(magrittr)
library(ggpubr)
library(jsonlite)

# import function
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# import validation peptide information
validation_info<-
  read_csv('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/sequence_validation.csv')

# import screen data
signif_data <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/signif_score_welch20_filter.csv')

# see if more double or single mutant lead to change in specificity

signif_data%>%
  dplyr::filter(condition == 'MTX' & side == 'stronger' & pool  %in% c('P4','P5', 'P6'))%>%
  dplyr::mutate(PBD=case_when(pool == 'P6' ~ '385',
                       pool == 'P5' ~ '366', 
                       pool == 'P4' ~ '363'))->specificity_data


specificity_data%<>%
  dplyr::filter(family != PBD)


specificity_data%<>%
  dplyr::mutate(mutant = case_when(max != deg ~ 'single', 
                            max1 != deg1 & max2==deg2 ~ 'single', 
                            max1 == deg1 & max2!=deg2 ~ 'single', 
                            max1 != deg1 & max2!=deg2 ~ 'double'))
specificity_data%>%
dplyr::mutate(PBD=case_when(pool == 'P6' ~ 'PBD 385',
                            pool == 'P5' ~ 'PBD 366', 
                            pool == 'P4' ~ 'PBD 363'))->specificity_data

ggplot(specificity_data)+
  facet_grid(cols = vars(PBD), scales = 'free_x')+
  geom_bar(aes(x = family, fill = mutant), color ='black', 
           position = position_dodge2(preserve = 'single'))+
  scale_fill_manual(values =c('grey30','grey80'))+
  ylab('count')+
  xlab('peptide family')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())->single_v_double

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/specificity_loss.png', 
       width = 5, height = 3)


saveRDS(single_v_double, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/specificity_loss_single_v_double.rds')


specificity_data%>%
  select(pos1, pos2, pos, aa_seq, PBD, med_norm, family, mutant)%>%
  filter(mutant == 'single')%>%
  pivot_longer(values_to = 'position', cols = c(pos1, pos2, pos), values_drop_na = TRUE)->specificity_position


ggplot(specificity_position)+
  facet_grid(cols = vars(PBD))+
  geom_bar(aes(x = position, fill = family), color ='black', linewidth = 0.4, width = 0.6, lineend = 'square')+
  scale_fill_manual(values =c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  ylab('count')+
  xlab('position')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank()) ->position_v_specificity

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/specificity_v_position.png', 
       width = 5, height = 3)

saveRDS(position_v_specificity, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/specificity_loss_v_position.rds')

# visualise examples
signif_data%>%
  filter(condition == 'MTX' & pool  %in% c('P4','P5', 'P6'))%>%
  mutate(PBD=case_when(pool == 'P6' ~ '385',
                       pool == 'P5' ~ '366', 
                       pool == 'P4' ~ '363'))->specificity_data

specificity_data%<>%
  mutate(mutant = case_when(max != deg ~ 'single', 
                            max1 != deg1 & max2==deg2 ~ 'single', 
                            max1 == deg1 & max2!=deg2 ~ 'single', 
                            max1 != deg1 & max2!=deg2 ~ 'double', 
                            max1 == deg1 & max2 == deg2 ~ 'PWMmax', 
                            max == deg ~ 'PWMmax'))


specificity_data%>%
  filter(mutant == 'single' & is.na(pos))%>%
  dplyr::mutate(pos = case_when(   max1 != deg1 & max2==deg2 ~ pos1, 
                                   max1 == deg1 & max2!=deg2 ~  pos2), 
                deg = case_when(   max1 != deg1 & max2==deg2 ~ deg1, 
                                   max1 == deg1 & max2!=deg2 ~  deg2), 
                max = case_when(   max1 != deg1 & max2==deg2 ~ max1, 
                                   max1 == deg1 & max2!=deg2 ~  max2))->single_double

specificity_data%<>%
  filter(!is.na(pos))%>%
  bind_rows(single_double)


specificity_data%>%
  dplyr::filter(mutant %in% c('single', 'PWMmax'))%>%
  select(aa_seq, deg, pos, max, PBD, med_norm, mutant, family)%>%
  mutate(mutation = paste0(max, pos, deg), 
         end_aa =  paste0(pos, deg), 
         start_aa =  paste0(max, pos))->mut_spe


mut_spe%>%
  filter(PBD == 363)%>%
  filter(end_aa %in% c('6D', '5W', '5W') | mutant =='PWMmax')->test


test%>%
  filter(mutation %in% c('F6D', 'H5W', 'W9W'))->sub_test

test%>%
  filter(end_aa == '6D' |  mutant =='PWMmax')%>%
ggplot()+
  geom_point(aes(mutant, med_norm, color = family, shape = mutant), size = 2.5)+
  geom_line(aes(mutant, med_norm, group =family, color =family), linetype = 'dashed')+
  geom_text(data = test[test$mutant == 'PWMmax' & test$pos == 6, ], 
            aes(mutant, med_norm+0.05, label = start_aa), size = 3)+
  scale_x_discrete(limits = c('PWMmax', 'single'), 
                   labels =c(expression(paste('PWM'[max])), 'X6D'))+
  scale_color_manual(values =c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  scale_shape_manual(values = c(16,4))+
  ylab('PPI score ~ PBD 363')+
  xlab('peptide')+
  ylim(-0.1,1)+
  t+
  theme(legend.position = 'none')->ex_1

test%<>%
  filter(mutant =='PWMmax')%>%
  mutate(mutant = 'single')%>%
  bind_rows(test)

test%>%
  filter(end_aa == '5W' |  mutant =='PWMmax')%>%
  ggplot()+
  geom_point(aes(mutant, med_norm, color = family, shape = mutant), size =2.5)+
  geom_line(aes(mutant, med_norm, group =family, color =family), linetype = 'dashed')+
  geom_text(data =  test[test$mutant == 'PWMmax' & test$pos == 5, ], 
            aes(mutant, med_norm+0.05, label = start_aa), size = 3)+
  scale_shape_manual(values = c(16,4))+
  scale_x_discrete(limits = c('PWMmax', 'single'), 
                   labels = c(expression(paste('PWM'[max])), 'X5W'))+
  scale_color_manual(values =c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  ylab('PPI score ~ PBD 363')+
  xlab('peptide')+
  ylim(-0.1,1)+
  t+
  theme(legend.position = 'none')->ex_2



exemples<-plot_grid(ex_2, ex_1)
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/exemples_spe.png', exemples)

saveRDS(ex_1, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ex1_spe.rds')
saveRDS(ex_2, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ex2_spe.rds')


# filter screen data to keep only the peptide present in the validation subset
# and in PCA condition

signif_data%>%
  filter(aa_seq %in% validation_info$aa_seq & condition == 'MTX')->data_validated


data_validated%<>%
  left_join(validation_info[, -1],
            join_by(aa_seq))

# select for PPI assay and specificity loss vs reference peptides

data_validated%>%
  filter(validation %in% c('specificity loss', 'reference') & pool %in% c('P4', 'P5', 'P6'))->specificity_data


specificity_data%>%
  filter(side == 'stronger' | c(validation == 'reference'))%>%
  select(aa_seq, med_norm, validation, family, pool)%>%
  unique()-> sub_vis

specificity_data%>%
  filter(side == 'stronger' | c(validation == 'reference'))%>%
  filter(aa_seq %in% unlist(sub_vis[sub_vis$validation == 'specificity loss', 'aa_seq']))%>%
  select(aa_seq, deg, pos, max, pool)%>%
  mutate(mutation = paste0(max, pos, deg))->mut_spe

sub_vis%<>%
  left_join(mut_spe[, c('aa_seq', 'mutation', 'pool')])

sub_vis$pool<-
  factor(sub_vis$pool,
         levels = c('P4', 'P5', 'P6'), 
         labels = c('PBD 363', 'PBD 366', 'PBD 385'))


sub_vis%>%
  filter(validation=='reference')%>%
select(family, med_norm, pool)->ref_pep

sub_vis%<>%
  filter(validation!='reference')%>%
  left_join(ref_pep, join_by(family, pool), suffix = c('.spe', '.ref'))

ref_pep%>%
  filter((pool == 'PBD 363' & family == 363) |
           (pool == 'PBD 366' & family == 366)|
           (pool == 'PBD 385' & family == 385))->sub_ref

sub_vis%>%
ggplot()+
  facet_wrap(vars(pool), scales ='free_y')+
  geom_linerange(aes(y= mutation, xmin = med_norm.ref, xmax =med_norm.spe), color = 'black')+
  geom_point(aes(y= mutation, med_norm.ref, color =family), shape =19)+
  geom_point(aes(y=mutation, med_norm.spe), shape =17, color = 'black')+
  geom_vline(data =sub_ref, 
             aes(xintercept = med_norm, color =family), linetype = 'dashed')+
  scale_color_manual(values =c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  xlim(0,1)+
  labs(x = 'PPI score', y = 'substitution vs PWMmax')+
  t











