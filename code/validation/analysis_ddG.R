# AUTHOR : Pascale Lemieux
# Comparison of ddG values obtained with FoldX and RosettaFlexddG with PCA scores

library(ggpubr)
library(tidyverse)
library(rstatix)
library(magrittr)


setwd('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/data/validation/')
source('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/code/functions.R')

all_ddg<-
  read_csv('all_ddg.csv')

exp_score<-
  read_csv('../PCA2/signif_score_welch20_filter.csv')

exp_score%<>%
  filter(condition=='MTX' & assay =='PPI')%>%
  filter( (pool =='P4' & family == '363')|
            (pool =='P5' & family == '366')|
            (pool =='P6' & family == '385'))

exp_score%>%
  #filter(method == 'foldx')%>%
  filter(side == 'stronger')%>%
  summarise(n())


full_join(all_ddg, exp_score[, c('aa_seq', 'med_norm', 'peptide_max', 'side', 'family')], 
          join_by(sequence == aa_seq))->comp_ddg_exp


comp_ddg_exp[comp_ddg_exp$peptide_max, 'ddg_binding']<-0


comp_ddg_exp%>%
  mutate(side = factor(comp_ddg_exp$side, 
                       levels = c('stronger', 'no significant difference', 'weaker')))%>%
  filter(method=='foldx')%>%
ggplot()+
  facet_grid(cols = vars(family.y))+
  geom_point(aes(ddg_binding, med_norm, color = side, shape = side, alpha=side))+
  stat_cor(aes(ddg_binding, med_norm), method = 'spearman', label.sep = '\n', 
           cor.coef.name = 'rho', p.digits = 3, label.x = 5,  size = 2.5)+
  scale_color_manual(values = c( '#4e67c8ff', 'grey', '#373737'))+
  scale_shape_manual(values = c(19,1,1))+
  scale_alpha_manual(values = c(1,0.3, 1))+
  geom_point(data = comp_ddg_exp[comp_ddg_exp$peptide_max, ], 
             aes(ddg_binding, med_norm), color ='black')+
  xlab(expression(paste(Delta, Delta, 'G binding (foldX)')))+
  ylab('PPI score')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())->Dpanel

comp_ddg_exp%>%
  mutate(side = factor(comp_ddg_exp$side, 
                       levels = c('stronger', 'no significant difference', 'weaker'), 
                       labels = c('stronger', 'no significant\ndifference', 'weaker')))%>%
  filter(method=='foldx')%>%
ggplot()+
  #geom_vline(aes(xintercept = 0), linetype ='dashed')+ 
  geom_violin(aes(x = side, y = ddg_binding, fill =side), 
              quantile.colour = 'black', quantile.linetype = 'solid', trim =F, alpha=0.5)+
  scale_fill_manual(values = c( '#4e67c8ff', 'grey', '#373737'))+
  ylab(expression(paste(Delta, Delta, 'G binding (foldX)')))+
  #stat_compare_means(aes(x = side, y = ddg_binding), method = 't.test', ref.group = 'no significant difference')+
  t+
  theme(legend.position = 'none', 
        axis.title.x = element_blank())->Cpanel


comp_ddg_exp%>%
  drop_na(method)%>%
  select(family.y, sequence, mutations, med_norm, ddg_binding, method, side)%>%
  pivot_wider(values_from = ddg_binding, names_from = method)->test

test%>%
  mutate(side = factor(test$side, 
                       levels = c('stronger', 'no significant difference', 'weaker')))%>%
ggplot()+
 
  geom_point(aes(foldx, flexddg), shape=1, alpha=0.4)+
  #scale_color_manual(values = c( '#4e67c8ff', 'grey', '#373737'))+
  #scale_shape_manual(values = c(19,1,1))+
  #scale_alpha_manual(values = c(1,0.4, 0.4))+
   geom_hline(aes(yintercept=1), color ='red', linetype='dashed')+
  geom_vline(aes(xintercept=1), color ='red', , linetype='dashed')+
  stat_cor(aes(foldx, flexddg), method = 'spearman', label.sep = '\n', 
           cor.coef.name = 'rho', size = 2.5)+
  xlab(expression(paste(Delta, Delta, 'G binding (foldX)')))+
  ylab(expression(paste(Delta, Delta, 'G binding (flexddg)')))+
  lims(x =c(-1.5, 12), y=c(-1.5,12))+
  t+
  theme(legend.position = 'none')->Bpanel


library(cowplot)

top<-plot_grid(Bpanel, Cpanel, align = 'h', axis = 'tb', labels = c('B', 'D'), 
          label_size = 12, label_fontface = 'plain')

plot_grid(top, Dpanel, nrow = 2, rel_heights = c(0.8,1),  labels = c('', 'C'), 
          label_size = 12, label_fontface = 'plain')



comp_ddg_exp%>%
  filter(ddg_binding<0 & (ddg_binding + ddg_binding_sd ) < 0)%>%
  filter(method=='foldx')->signif_stab

signif_stab%>%
  ggplot()+
  geom_point(aes(y=med_norm, side))


table(signif_stab$side)
16/(16+17+167)

table(comp_ddg_exp$side)
124/(2045+124+1281)


comp_ddg_exp%>%
  filter(method == 'foldx')%>%
  filter(side == 'stronger')%>%
  summarise(n())

# PPI stronger in silico (foldx) and in vivo
a = comp_ddg_exp%>%
  filter(method == 'foldx')%>%
  filter(side == 'stronger'& ddg_binding+ddg_binding_sd < 0)%>%
  summarise(n())%>%
  pull('n()')

# PPI stronger in silico but not in vivo
b = comp_ddg_exp%>%
  filter(method == 'foldx')%>%
  filter(side != 'stronger' & ddg_binding+ddg_binding_sd < 0)%>%
  summarise(n())%>%
  pull('n()')

# PPI stronger in vivo but not in silico
c = comp_ddg_exp%>%
  filter(method == 'foldx')%>%
  filter(side == 'stronger'& ddg_binding+ddg_binding_sd > 0)%>%
  summarise(n())%>%
  pull('n()')
  
# PPI not stronger in vivo and in silico
d = comp_ddg_exp%>%
  filter(method == 'foldx')%>%
  filter(side != 'stronger'& ddg_binding+ddg_binding_sd > 0)%>%
  summarise(n())%>%
  pull('n()')

f_matrix<-matrix(nrow = 2, ncol = 2, 
                 data = c(a,b,c,d), byrow = T)  
  
fisher_test(xtab = f_matrix, alternative = 'greater')
# p-value = 0.00258, signif association

opt_ddg<-
  read_csv('opt_ddg.csv')

gc_opt<-
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/data/validation/gc_score.csv')


opt_ddg%<>%
  left_join(gc_opt, join_by(sequence==peptide))

opt_ddg%>%
  filter(method=='foldx')%>%
ggplot()+
  geom_boxplot(aes(x =as.factor(family), y = ddg_binding, fill = purpose),
              position = position_dodge2(width=0.2, padding = 0.2), alpha = 0.6)+
  ylab(expression(paste(Delta, Delta, 'G binding (foldX)')))+
  xlab('peptide family')+
  scale_fill_manual(values = c('#4e67c8ff','#A06AB4'), 
                     labels = c('stronger\npeptides', 'optimized\npeptides'))+
  t+
  theme(legend.position = 'top', 
        legend.title = element_blank())->Epanel


opt_ddg%>%
  filter(method=='foldx')%>%
  ggplot()+
  geom_point(aes(x =ppi_exp, y = ddg_binding, color = as.factor(family), shape =purpose),
               position = position_dodge2(width=0.2, padding = 0.2))+
  stat_cor(aes(x =ppi_exp, y = ddg_binding), method ='pearson')+
  ylab(expression(paste(Delta, Delta, 'G binding (foldX)')))+
  xlab('peptide family')+
  scale_color_manual(values = c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  #scale_fill_manual(values = c('#4e67c8ff','#A06AB4'), 
   #                 labels = c('stronger\npeptides', 'optimized\npeptides'))+
  t+
  theme(legend.position = 'top', 
        legend.title = element_blank())#->Epanel


top<-plot_grid(Bpanel, Cpanel, Epanel, align = 'h', axis = 'tb', labels = c('B', 'D', 'E'), 
               label_size = 12, label_fontface = 'plain', nrow =1)

plot_grid(top, Dpanel, nrow = 2, rel_heights = c(0.8,1),  labels = c('', 'C'), 
          label_size = 12, label_fontface = 'plain')

library(svglite)
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/in_silico.svg', 
       width = 8, height = 5)


