# AUTHOR : Pascale Lemieux
# Comparison of bulk competition screen results with individual validation
library(tidyverse)
library(ggpubr)
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# validation avail
summary_avail<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/validation/condensed_validation_nov2025.csv')
pep_reference<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/pep_reference.csv')
pep_validation<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/peptide_validation.csv')

summary_avail%<>%
  data_summary( 'auc', 'D3')


# import median scores of the July2025 experiment (single mutant availability)
specificity <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/filtered_coefficient.csv')

avail_screen2 <- subset(specificity, pool == 'P1' & condition == 'MTX')

# import the scores from the February2025 experiment (double mutant availability)
avail <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/norm_avail_sel_coef.csv')


##here combine both screens to compare with the validation growthcurves
avail%>%
  filter(assay == 'Avail')%>%
  select(family, aa_seq,norm_sel_coef)%>%
  data_summary('norm_sel_coef', 'aa_seq')->avail_1

avail_screen2%>%
  filter(condition =='MTX')%>%
  select(aa_seq, norm_sel_coef, family)%>%
  data_summary('norm_sel_coef', 'aa_seq')->avail_2


bind_rows(avail_1, avail_2, .id = 'screen')->all_avail

summary_avail%>%
  left_join(pep_reference, join_by(D3==ID))->summary_avail

summary_avail[!is.na(summary_avail$`aa sequence`), 'D3']<- summary_avail[!is.na(summary_avail$`aa sequence`), 'aa sequence']

right_join(all_avail, summary_avail, join_by(aa_seq == D3), suffix = c('.score', '.auc'))->summary_avail

summary_avail%<>%
  mutate(purpose = case_when(aa_seq %in% pep_reference$`aa sequence` ~ 'reference', 
                            aa_seq %in% pep_validation$aa_seq ~ 'validation'))


# comparison availability gc and screen
summary_avail%>%
  drop_na(screen)%>%
  ggplot(aes(auc, norm_sel_coef))+
  #facet_grid(cols =vars(screen))+
  geom_pointrange(aes(xmin =auc-sd.auc ,xmax = auc+sd.auc), color = 'darkgrey', size = 0.1, linewidth = 0.8, alpha= 0.5)+ 
  geom_pointrange(aes(ymin = norm_sel_coef-sd.score, ymax = norm_sel_coef+sd.score), 
                  color = 'darkgrey', size = 0.1, linewidth = 0.8, alpha= 0.5)+
  geom_point(aes(color =purpose), size = 2)+
  ylim(-0.2, 1.2)+
  #geom_text(aes(label = D3))+
  stat_cor(method = 'spearman', label.sep = '\n', cor.coef.name = 'rho', size = 3)+
  #scale_color_viridis_c()+
  scale_color_manual(values =c('red','grey'))+
  xlab('avail. signal \n (corrected AUC)')+
  ylab('avail. score')+
  t+
  theme(legend.position = 'bottom',
        legend.title = element_blank())->val_avail

saveRDS(val_avail, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/validation_avail.rds')

PPI_val<-
read_csv( '~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/val_ppi.csv')

# validation of increased PPI for max PPI peptide 

PPI_val%>%
  filter(validation == 'max PPI' | method == 'PWM_top')%>%
  filter(!(ID %in% c('sh3.1', 'pdz.1')))%>%
  filter(D3 != 'empty')%>%
  ggplot()+
  facet_grid(cols = vars(D12), scales = 'free')+
  geom_errorbar(aes(x = D3, y = auc, ymax = auc+sd.auc, ymin = auc-sd.auc))+
  geom_point(aes(x = D3, y = auc, color = factor(family), shape = method), size = 2)+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF")[5:7])+
  scale_shape_manual(values = c(17,16))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle =90, hjust =1))+
  guides(color = guide_legend(title = 'pep. family'))-> max_val
saveRDS(max_val, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/validation_max.rds')

ggsave('max_peptide_validation.png')

PPI_val%>%
  filter(validation == 'specificity loss' | is.na(validation))%>%
  filter(!(ID %in% c('sh3.1', 'pdz.1')))%>%
  filter(D3 != 'empty')%>%
  drop_na(family)%>%
  ggplot()+
  facet_grid(cols = vars(D12), scales = 'free')+
  geom_errorbar(aes(x = D3, y = auc, ymax = auc+sd.auc, ymin = auc-sd.auc))+
  geom_point(aes(x = D3, y = auc, color = factor(family), shape = method), size =2)+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF")[5:7])+
  scale_shape_manual(values = c(17,16))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle =90, hjust =1))+
  guides(color = guide_legend(title = 'pep. family'))->spe_val

ggsave('specificity_peptide_validation.png')

saveRDS(spe_val, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/validation_specificity.rds')


PPI_val%>%
  filter(validation == 'PCA optimization' | is.na(validation))%>%
  #filter(!(ID %in% c('sh3.1', 'pdz.1')))%>%
  filter(D3 != 'empty')%>%
  ggplot()+
  facet_grid(cols = vars(D12), scales = 'free')+
  geom_errorbar(aes(x = D3, y = auc, ymax = auc+sd.auc, ymin = auc-sd.auc))+
  geom_point(aes(x = D3, y = auc, color = factor(family), shape = method))+
  #scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF")[5:7])+
  scale_shape_manual(values = c(17,16))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle =90, hjust =1))+
  guides(color = guide_legend(title = 'pep. family'))


specificity%>%
  filter(pool %in% c('P4', 'P5', 'P6') & condition=='MTX')%>%
  select(aa_seq, norm_sel_coef, pool)%>%
  dplyr::mutate(D12 = case_when(pool == 'P4' ~ 363, 
                                pool == 'P5' ~ 366, 
                                pool == 'P6' ~ 385))->simple_spe
summary_spe<-tibble()

for (i in c(363, 366, 385)) {

    simple_spe%>%
    dplyr::filter(D12 == i)%>%
    data_summary('norm_sel_coef', 'aa_seq')->x
    x$D12<-i
    summary_spe<-bind_rows(summary_spe, x)
}


summary_spe%>%
  dplyr::mutate(D12 = as.character(D12))%>%
  right_join(PPI_val, join_by(D12 == D12, aa_seq==D3))->comp_screen

summary_spe[summary_spe$aa_seq %in% unlist(pep_reference[is.na(pep_reference$method), 'aa sequence']), ]->ref_pep

ref_pep%>%
  dplyr::mutate(D12 = as.character(D12))%>%
  right_join(PPI_val[, -1], join_by(aa_seq==D3))->comp_screen_ref


bind_rows(comp_screen, comp_screen_ref)->comp_screen
  


comp_screen%>%
  filter(strain == 'BY4741' & D12 %in% c('363', '366', '385'))%>%
  dplyr::mutate(D12 = str_c('PBD ', D12))%>%
  drop_na(family)%>%
  ggplot(aes(auc, norm_sel_coef))+
  facet_wrap(vars(D12))+
  geom_pointrange(aes(xmin = auc-sd.auc,xmax = auc+sd.auc), color = 'darkgrey', size = 0.1, linewidth = 0.8, alpha= 0.5)+ 
  geom_pointrange(aes(ymin = norm_sel_coef-sd, ymax = norm_sel_coef+sd), 
                  color = 'darkgrey', size = 0.1, linewidth = 0.8, alpha= 0.5)+
  geom_point(aes(color=family, shape = method))+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.sep = '\n')+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF",'red', 'red')[5:9],
                    labels = c('363', '366', '385','reference' ,'reference'))+
  scale_shape_manual(values = c(17,16), labels = c(expression(paste('PWM'[max])), 
                                                   'validation'))+
  xlab('PPI signal\n(mean corrected AUC)')+
  ylab('PPI score')+
  t+
  guides(color=guide_legend(title = 'peptide\nfamily'), 
         shape=guide_legend(title =''))#->val_screen

saveRDS(val_screen, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/validation_screen.rds')



comp_screen%>%
  select(aa_seq, D12, validation)%>%
  filter(validation !='PCA optimization')%>%
  unique()%>%
  nrow()
#214 PBD-peptide pairs


