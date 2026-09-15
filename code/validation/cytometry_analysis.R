# AUTHOR : Pascale Lemieux
# comparison of cytometry data with availability 

library(tidyverse)
library(ggpubr)
library(magrittr)

setwd(dir = '~/PL_projects/PL_papers/PPI_optimization_paper/data/')

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

val_pep<-
  read_csv('peptide_validation.csv')

# import cytometry data from may 2026
simple_GFP<-read_csv('validation/cytometry_data.csv')
colnames(simple_GFP)[3]<-'sd.med_GFP'

simple_GFP$seq<-
gsub('pGD110', 'negative', simple_GFP$seq)

simple_GFP$seq<-
  gsub('pPL9', 'empty', simple_GFP$seq)



#import and format avail validation
val_avail<-read_csv('~/PL_projects/growthcurves/validation_Nov2025/condensed_validation_nov2025.csv')

sum_avail<-
  data_summary(val_avail, 'auc', 'D3')

sum_avail$D3<-
gsub('pGD110', 'empty', sum_avail$D3)



# add PCA optimized avail data from april 2026

avail_opt<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/validation/condensed_avail_optimized_peptide.csv')

avail_opt%>%
  dplyr::filter(is.na(`note sequencing`))%>%
  data_summary('auc', 'peptide')->sum_avail_opt

colnames(sum_avail_opt)[1]<-'D3'

sum_avail<-bind_rows(sum_avail, sum_avail_opt, .id = 'batch')



# keep avail batch 1 for the reference sequences

sum_avail<-sum_avail[!duplicated(sum_avail$D3), ]


# merge avail and gfp

full_join(simple_GFP, sum_avail, 
          join_by(seq == D3))->avail_gfp

full_join(avail_gfp, val_pep, 
          join_by(seq == aa_seq))->avail_gfp

avail_gfp$validation<-
replace_na(avail_gfp$validation, 'reference')

avail_gfp$seq<-
  gsub('empty', 'linker', avail_gfp$seq)

library(ggrepel)


avail_gfp%>%
  dplyr::filter(!validation %in% c('stop codon', 'PCA optimization', 'specificity loss'))%>%
ggplot()+
  geom_errorbar(aes(x=auc, y=med_GFP, xmin =auc-sd, xmax = auc+sd), color = 'grey50', alpha = 0.6)+
  geom_errorbar(aes(x=auc, y=med_GFP, ymin =med_GFP-sd.med_GFP, ymax = med_GFP+sd.med_GFP), color = 'grey50', alpha = 0.6)+
  geom_point(aes(x=auc, y=med_GFP, color = validation))+
  stat_cor(aes(x=auc, y=med_GFP,  color = validation), size = 3, method = 'spearman', cor.coef.name = 'rho', show.legend = FALSE)+
  geom_text_repel(data = avail_gfp[avail_gfp$seq %in% c('299-1', '246-1', '363-1', '366-1', '385-1', 'linker', '250-1'), ],
                  aes(x=auc, y=med_GFP, label = seq),
                  min.segment.length = 0.01,
                  size =3, 
                  nudge_x = 0.3,
                  box.padding = 0.3,
                  nudge_y = -1,
                  segment.curvature = -0.6,
                  segment.ncp = 0.3,
                  segment.angle = 30
  )+
  scale_color_manual(values= c('#C3D37A','#373737'))+
  scale_x_continuous(transform = 'log2')+
  ylab('mean GFP-peptide \nfluorescence (a.u.)')+
  xlab('Availability \n(corrected AUC)')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())-> avail_vs_GFP

saveRDS(avail_vs_GFP, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/gfp_vs_avail.rds')
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/gfp_vs_avail.png', avail_vs_GFP)

avail_gfp%>%
  dplyr::filter(validation %in% c('stop codon'))%>%
  ggplot()+
  geom_errorbar(aes(x=auc, y=med_GFP, xmin =auc-sd, xmax = auc+sd), color = 'grey50', alpha = 0.6)+
  geom_errorbar(aes(x=auc, y=med_GFP, ymin =med_GFP-sd.med_GFP, ymax = med_GFP+sd.med_GFP), color = 'grey50', alpha = 0.6)+
  geom_point(aes(x=auc, y=med_GFP, color = validation))+
  stat_cor(aes(x=auc, y=med_GFP,  color = validation), size = 3, method = 'spearman', cor.coef.name = 'rho', show.legend = FALSE)+
  geom_text_repel(
                  aes(x=auc, y=med_GFP, label = seq),
                  min.segment.length = 0.01,
                  size =3, 
                  nudge_x = 0.3,
                  box.padding = 0.3,
                  nudge_y = -1,
                  segment.curvature = -0.6,
                  segment.ncp = 0.3,
                  segment.angle = 30
  )+
  #scale_color_manual(values= c('#C3D37A','#373737'))+
  scale_x_continuous(transform = 'log2')+
  ylab('mean GFP-peptide \nfluorescence (a.u.)')+
  xlab('Availability \n(corrected AUC)')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())


avail_gfp%>%
  dplyr::filter(!validation %in% c('stop codon', 'specificity loss', 'reference') |
                  seq %in% c('363-1', '366-1', '385-1'))%>%
ggplot()+
  geom_boxplot(aes(x = validation, y = med_GFP), width = 0.2)+
  geom_jitter(aes(x = validation, y = med_GFP, color = validation), shape = 1, width = 0.1)+
  scale_color_manual(values = c('#C3D37A' , '#e3b8ffff','#8FFFD9'))+
  #stat_compare_means(aes(x = validation, y = med_GFP, group = validation), na.rm = T, 
   #                  method = 'wilcox.test', ref.group = 'reference')+
  ylab('mean GFP-peptide \nfluorescence (a.u.)')+
  xlab('peptide type')+
  t+
  theme(legend.position = 'none')->PCA_opt_GFP


saveRDS(PCA_opt_GFP, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/PCAopt_GFP.rds')
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/PCAopt_GFP.png', PCA_opt_GFP)

avail_gfp%>%
  dplyr::filter(!validation %in% c('stop codon', 'specificity loss', 'reference') |
                  seq %in% c('363-1', '366-1', '385-1'))%>%
  ggplot()+
  geom_boxplot(aes(x = validation, y = auc), width = 0.2)+
  geom_jitter(aes(x = validation, y = auc, color = validation), shape = 1, width = 0.1)+
  scale_color_manual(values = c('#C3D37A' , '#e3b8ffff','#8FFFD9'))+
  stat_compare_means(aes(x = validation, y = auc, group = validation), na.rm = T, 
                     method = 'wilcox.test', ref.group = 'reference')+
  ylab('mean GFP-peptide \nfluorescence (a.u.)')+
  xlab('peptide type')+
  #scale_x_discrete(limits = c('reference', 'max PPI', 'PCA optimization'), 
  #                 labels = c((expression(paste('PWM'[max]))), 'stronger\nvalidation', 'opt. peptides'))+
  t+
  theme(legend.position = 'none')


avail_gfp[avail_gfp$seq == '363-1', 'family']<-363
avail_gfp[avail_gfp$seq == '366-1', 'family']<-366
avail_gfp[avail_gfp$seq == '385-1', 'family']<-385



avail_gfp%>%
  dplyr::filter(validation %in% c('reference', 'max PPI', 'PCA optimization') & !is.na(family))%>%
  dplyr::group_by(family, validation)%>%
  dplyr::summarise(med_auc_family = median(auc, na.rm = TRUE), 
                   med_GFP_family = median(med_GFP, na.rm = TRUE))%>%
  pivot_wider(values_from = c(med_auc_family, med_GFP_family), names_from = c(validation))->median_abundance

median_abundance%>%
  dplyr::mutate(family = as.factor(family))%>%
ggplot()+ 
  geom_segment(aes(y=family, x = `med_GFP_family_reference`, xend = `med_GFP_family_max PPI`))+
  geom_point(aes(y=family, x = `med_GFP_family_max PPI`),shape =4, size = 3)+
  geom_point(aes(y=family, x = `med_GFP_family_reference`, color = family), size =3)+
  scale_color_manual(values = c( "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  xlab('median GFP signal\n(a.u.)')+
  ylab('peptide family')+
  t+
  theme(legend.position = 'none')->GFP_median_grah


median_abundance%>%
  dplyr::mutate(family = as.factor(family))%>%
  ggplot()+ 
  geom_segment(aes(y=family, x = `med_auc_family_reference`, xend = `med_auc_family_max PPI`))+
  geom_point(aes(y=family, x = `med_auc_family_max PPI`), shape = 4, size = 3)+
  geom_point(aes(y=family, x = `med_auc_family_reference`,  color = family), size =3)+
  scale_color_manual(values = c( "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  xlab('median avail. signal\n(corrected AUC)')+
  ylab('peptide family')+
  t+
  theme(legend.position = 'none')->Avail_median_grah

saveRDS(GFP_median_grah, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/med_GFP_graph.rds')
saveRDS(Avail_median_grah, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/med_Avail_graph.rds')
