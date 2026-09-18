# AUTHOR : Pascale Lemieux
# Figures assembly

library(tidyverse)
library(viridisLite)
library(cowplot)
library(svglite)
library(magick)
library(ggpubr)

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/')
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# Figure 1

Fig2A <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2A.rds')
Fig2B <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2B.rds')

Fig2C <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2C.rds')+
  scale_color_manual(values = c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))

Fig2D <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/fig2D.rds')+
  scale_fill_manual(values = c('#373737', 'grey', '#4e67c8ff' ))



Fig2E <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig2E.rds')

Fig2F <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2F.rds')

Fig2G <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2G.rds')+
    scale_color_manual(values = c( 'grey', '#4e67c8ff','#373737' ))
  

 p1 <- 
   ggdraw()+
   draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig1.png')

 panelAB<-
  plot_grid(Fig2A, 
            Fig2B, ncol= 2, align = 'h', axis = 'tb', labels = c('D', 'E'), 
            label_fontface = 'plain', label_size = 14)

panelDF<-
 plot_grid(panelAB, Fig2C, nrow =2,  labels = c('', 'F'), 
           label_fontface = 'plain', label_size = 14, rel_heights = c(0.5, 1))

panelAD <- 
  plot_grid(p1, panelDF,ncol = 2,
            rel_widths = c(1.2, 1), labels = c('', ''), label_fontface = 'plain', label_size = 14)

PPI_leg<-get_legend(Fig2E)
binding_leg<-get_legend(Fig2G)

leg<-plot_grid( PPI_leg, binding_leg, rel_widths = c(2,1))



panelEG <- 
  plot_grid(
    Fig2E+theme(legend.position = 'none'), 
    Fig2F+theme(legend.position = 'none'),
    Fig2G+theme(legend.position = 'none'),
    nrow = 1, ncol = 3, axis = 'tb', align = 'h',
    labels = c('G', 'H', 'I'), rel_widths = c(1,1,1), label_fontface = 'plain', label_size = 14)


Fig1 <- 
  plot_grid(panelAD, panelEG, leg,
            nrow = 3, rel_heights = c(1.2, 1, 0.1))


#ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig1_v2.svg', height = 7, width = 10)


# Figure 2

Fig3A <- 
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig3A.png')

Fig3B<-
ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig3B.png')


Fig3C <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/validation_btw_screen.rds')

panelAC<-
plot_grid(Fig3A, Fig3B, Fig3C, rel_widths = c(1.4,0.85,0.71), ncol = 3, 
          labels = c('A', 'B', 'C'), label_fontface = 'plain', label_size = 14)

P4P5<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4P5.rds')+
  ylab('PPI score ~ PBD 366')+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position = 'bottom')+
  geom_vline(aes(xintercept = 0.25), linetype ='dashed')+
  geom_hline(aes(yintercept = 0.25), linetype ='dashed')+
  scale_alpha_manual(values = c(0,1), labels = c('', 'validation'))+
  guides(color = guide_legend(title = 'peptide family'), 
         alpha = guide_legend(title = ''))

leg_comp_spe<-get_legend(P4P5)


P4P5<-
P4P5+
  theme(legend.position = 'none')+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
                     limits = c(1,-0.1))+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
                     limits = c(1,-0.1))

P4P6<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4P6.rds')+
  ylab('PPI score ~ PBD 385')+
  xlab('PPI score ~ PBD 363')+
  geom_vline(aes(xintercept = 0.25), linetype ='dashed')+
  geom_hline(aes(yintercept = 0.25), linetype ='dashed')+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
                     limits = c(1,-0.1))+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
                     limits = c(1,-0.1))
  
P5P6<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P5P6.rds')+
  xlab('PPI score ~ PBD 366')+
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())+
  geom_vline(aes(xintercept = 0.25), linetype ='dashed')+
  geom_hline(aes(yintercept = 0.25), linetype ='dashed')+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
                     limits = c(1,-0.1))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
                     limits = c(1,-0.1))



Fig3F<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/pos_prop.rds')

Fig3E <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/comp_stronger_pep.rds')

  
Fig3H <- 
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ppi_all_vs_PWM.rds')


Fig3I1<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ex1_spe.rds')
Fig3I2<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ex2_spe.rds')


l_pos<-get_legend(Fig3E)

Fig5A <- 
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig5A.png')

pred_opt<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/ex_PCA_opt.rds' )+
  scale_color_manual( values = c('#4e67c8ff','#A06AB4','#373737'))

l_opt<-get_legend(pred_opt+theme(legend.position = 'bottom'))

Fig3I<-
plot_grid(Fig3I2, Fig3I1+theme(axis.text.y =element_blank(),
                                       axis.ticks.y = element_blank(),
                                       axis.title.y = element_blank()), 
          rel_widths = c(1, 0.85),
          labels = c( '', ''),  
          label_fontface = 'plain', label_size = 14)

Fig3I_l<-plot_grid(Fig3I, l_opt, nrow = 2, rel_heights = c(1, 0.3))

Fig3FG<-
plot_grid(Fig3F+theme(legend.position = 'bottom'),
         Fig3H, 
          align = 'h', axis = 'tb',
          nrow = 1, rel_widths = c(1,1),  labels = c('F', 'G'), 
          label_fontface = 'plain', label_size = 14)

Fig3FGH<-
plot_grid(Fig3FG, Fig3I_l, rel_widths = c(1.3, 0.93))

Fig5IJ<-
plot_grid(Fig5A, pred_opt+theme(legend.position = 'none'), 
          nrow = 1, rel_widths = c(1,1.2),  labels = c('I', 'J'), 
          label_fontface = 'plain', label_size = 14)



plot_grid(Fig3E+theme(legend.position = 'none'), 
Fig3FGH,
labels = c('E', ''), label_fontface = 'plain', label_size = 14,
nrow = 1,
rel_widths = c(0.7,2.3))->EH



Fig3DI<-
plot_grid(P4P5,
          EH,
          labels = c('D', ''), label_fontface = 'plain', label_size = 14,
          nrow = 1,
          rel_widths = c(1.1,3.3))


Fig3Db<-
plot_grid(P4P6, P5P6, rel_widths = c(1,0.85))



Fig3DH<-
plot_grid(Fig3Db, Fig5IJ, rel_widths = c(2.5, 2.9))


panelDI<-
plot_grid(Fig3DI, Fig3DH, nrow = 2, rel_heights = c(0.9, 1))

top_r<-
plot_grid( Fig5A, pred_opt+theme(legend.position = 'none'),
           nrow = 2,  labels = c( 'C' , 'D'), label_fontface = 'plain', label_size = 14)


top<-
  plot_grid(Fig3A, Fig3F+theme(legend.position = 'top'), 
            Fig5A, pred_opt+theme(legend.position = 'top'),
            rel_widths = c(0.8, 0.5, 0.5, 0.7), ncol = 4, 
            labels = c('A', 'B', 'C', 'D'), label_fontface = 'plain', label_size = 14)


middle<-
 plot_grid(Fig3B, Fig3C,  P4P5, Fig3E+theme(legend.position = 'none'), nrow = 1, 
           rel_widths = c(0.8, 0.6, 0.6, 0.5), ncol = 4, 
           labels = c('E', 'F', 'G', 'H'), label_fontface = 'plain', label_size = 14) 


bottom<-
  plot_grid(Fig3H+theme(legend.position = 'top'), 
            Fig3I, 
            P4P6, P5P6, ncol = 4, rel_widths = c(0.7, 0.7, 0.6, 0.5), 
            labels = c('I', 'J', '', ''), label_fontface = 'plain', label_size = 14)



Fig2<-
  plot_grid(top, middle, bottom, nrow = 3, rel_heights = c(0.9,1,1)) 

#ggsave(Fig3, filename ='~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig2v2.svg',
 #     width = 14, height =8.5)


### Fig 4 : availability

#FigA : availability schematics
Fig4A <- 
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig4A_1.svg')

# check quality data
Fig3B<-readRDS(file = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig3B.rds')

Fig3A<-readRDS(file = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig3A.rds')
Fig3A<- grid.grabExpr(grid.draw(Fig3A))

rep_check<-
plot_grid(Fig3A, Fig3B, ncol = 2, labels = c('B', 'C'), nrow =1,
          label_fontface = 'plain', label_size = 14)


Fig4ABC<-plot_grid(Fig4A, rep_check, labels = c('A', ''), nrow =2, rel_heights = c(1,0.7),
                   label_fontface = 'plain', label_size = 14)

# Fig4B : availability distribution (showing that PWM max is worse than average)
Fig4B<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig_avail_dist.rds')+
  scale_color_manual(values = c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))

# Fig4CD : Ppi avail association
stronger_ppi<-
  readRDS( '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/stronger_ppi_vs_avail.rds')+
    scale_fill_manual(values = c( '#373737', 'grey', '#4e67c8ff' ))+
  t+
  ggtitle('stronger\npeptides')+
  theme(legend.position = 'bottom', plot.title = element_text(hjust=0.5, size = 10), 
        legend.title = element_blank(), panel.border = element_blank(), 
        axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())


leg_avail_ppi<-
  cowplot::get_legend(stronger_ppi)

stronger_ppi<-
stronger_ppi+
  xlab('peptide family')+
  theme(legend.position = 'none')
 


all_ppi<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/all_ppi_vs_avail.rds')+
  #ggtitle('all PPIs')+
  xlab('peptide family')+
  ggtitle('all\npeptides')+
  scale_fill_manual(values = c( '#373737', 'grey', '#4e67c8ff' ))+
  theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = 'none', 
        panel.border = element_blank(), 
        axis.ticks.x = element_blank())

Fig4EF<-plot_grid(all_ppi, stronger_ppi, labels = c('E', 'F'), ncol =2,
                  label_fontface = 'plain', label_size = 14, rel_widths = c(1, 0.7))

FigEF_l<-plot_grid(Fig4EF, leg_avail_ppi, nrow =2, rel_heights = c(1, 0.2))

# Fig4F avail vs GFP-abundance
Fig4G<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/gfp_vs_avail.rds')+
  scale_color_manual(values = c('#4e67c8ff','#373737'), 
                     labels = c( 'stronger peptides','references'))+
  xlab('availability signal\n(corrected AUC)')+
  theme(legend.position = 'top')


stop_effect<-
readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/stop_avail_ppi.rds')+
  theme(legend.position = 'top')+
  xlab('PPI score \n  ')

GFP_med<-
readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/med_GFP_graph.rds')+
  xlab('median GFP-peptide\nsignal (a.u.)')

Avail_med<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/med_Avail_graph.rds')

opt_GFP<-
    readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/PCAopt_GFP.rds')+
    scale_color_manual( values = c('#4e67c8ff','#A06AB4','#373737'))+
    scale_x_discrete(limits = c('reference', 'max PPI', 'PCA optimization'), 
                     labels =  c((expression(paste('PWM'[max]))),'stronger\npeptides', 'opt.\npeptides'))+
    theme(axis.text.x = element_text(angle = 0, hjust =0.5, vjust = -1.5))
  
opt_avail<-
    readRDS( '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/PCAopt_Avail.rds')+
    scale_color_manual( values = c('#4e67c8ff','#A06AB4','#373737'))+
   scale_x_discrete(limits = c('reference', 'max PPI', 'PCA opt'), 
                   labels =  c((expression(paste('PWM'[max]))),'stronger\npeptides', 'opt.\npeptides'))+
  theme(axis.text.x = element_text(angle = 0, hjust =0.5, vjust = -1.5))

sticky_vs_avail<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/avail_vs_hydro.rds')+
  xlab('hydrophobicity score\n(cornell scale)')+
  ylab('availability\nscore')+
  theme(legend.position = 'none')

Fig4HK<-
  plot_grid(Avail_med, opt_avail, GFP_med, opt_GFP, 
            nrow = 2, ncol = 2, labels = c( 'H', 'J', 'I', 'K'), rel_widths = c(1, 1,1, 1), rel_heights = c(1,1,1,1),
            label_fontface = 'plain', label_size = 14, align = 'hv', axis = 'tblr')



Fig4ABCD<-
  plot_grid(Fig4ABC, Fig4B, FigEF_l, labels = c( '','D', '' ), 
            rel_widths = c(1.2, 1, 0.8), ncol =3,
            label_fontface = 'plain', label_size = 14)


FigEH<-
plot_grid(Fig4G, Fig4HK, sticky_vs_avail, ncol=3,
          rel_widths = c(1, 1.8, 1), label_fontface = 'plain', label_size = 14, labels = c('G', '', 'L'))

Fig4<-
plot_grid(Fig4ABCD, FigEH, nrow = 2, rel_heights = c(1, 1.2))

#ggsave(filename = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig4v2.svg', Fig3, height = 6, width=10)

# Figure S availability

SB<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/med_stop_position_effect_avail.rds')+
  ylab('availability score')
SA<-readRDS( '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/position_effect_avail.rds')+
  guides(fill =guide_colorbar(title='availability\nscore'))
SC<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/avail_vs_hydro_all.rds')

topS<-
plot_grid(SA, SB, 
          rel_widths = c(1,0.35), 
          align = 'h', axis = 'tb', 
          labels = c('A', 'B'), label_fontface = 'plain', label_size = 14)
S_avail<-
plot_grid(topS, SC, nrow = 2, rel_heights = c(1, 0.6), labels = c('', 'C'), 
          label_fontface = 'plain', label_size = 14)

ggsave(S_avail, filename ='~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS17.svg', 
       width = 10, height = 6)



# Figure S : processing PCA2#

T0_cov<-
  readRDS( '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/coverage_per_sample_T0_both.rds')+
  scale_x_continuous(breaks = c(1,10,100, 1000), transform = 'log10')
read_count<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/read_count_per_sample_family.rds')
sel_coef_timepoint<-
readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/sel_coef_distribution.rds')

norm_coef_distribution<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/norm_coef_distribution.rds')

norm_coef_distribution_bf<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/norm_coef_distribution_bf.rds')

tol_P5<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P5_tolerance.rds')+
  scale_color_manual(values = c( 'darkorange','#FFCC65','grey'))+
  guides(color = guide_legend(title = 'replication\nindex', position = 'bottom', ))



plot_grid(norm_coef_distribution_bf+
            theme(legend.position = 'none', axis.title.x = element_blank(), 
                  axis.text.x = element_blank(), 
                  axis.ticks.x =element_blank()), 
          norm_coef_distribution+
            theme(strip.background = element_blank(), 
                  strip.text = element_blank(), legend.position = 'none'), nrow =2, rel_heights = c(1, 1.2))->norm_sel_filt

plot_grid(T0_cov, tol_P5,  align = 'h', axis = 'tb',  labels = c('A', 'C'), 
          label_fontface = 'plain', label_size = 14, rel_widths = c(1, 0.7))->top_FS

plot_grid(top_FS, sel_coef_timepoint, norm_sel_filt, nrow = 3, align = 'v', axis = 'r',
          labels = c('', 'B', 'D'), 
          label_fontface = 'plain', label_size = 14)->FigS_PCA2_process

ggsave(FigS_PCA2_process, filename ='~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS5.svg', 
       width = 9, height =10)


# Figure S : replicability before vs after filtering

P1_filt<-
ggdraw()+
cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P1_norm_repcheck_filter.png')

P1<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P1_norm_sel_repcheck.png')


P4_filt<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P4_norm_repcheck_filter.png')

P4<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P4_MTX_norm_sel_repcheck.png')


P5_filt<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P5_norm_repcheck_filter.png')


P5<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P5_MTX_norm_sel_repcheck.png')


P6<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P6_MTX_norm_sel_repcheck.png')


P6_filt<-
  ggdraw()+
  cowplot::draw_image('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/P6_norm_repcheck_filter.png')


plot_grid(P1, P1_filt, 
          P4, P4_filt, 
          P5, P5_filt, 
          P6, P6_filt, 
          ncol = 2, nrow = 4)-> filter_MTX


ggsave(filter_MTX, filename = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS6.svg', 
       height = 12, width = 6)

# supp figure scores 
signif_count<-
readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/count_signif.rds')
stat_comp<-
readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/stat_sel_coef_20_filter.rds')+
  xlab('peptide family')+
  scale_alpha_manual(values = c(0.3, 1,1))

PWM_comp<-
readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/PWM_max_vs_variant.rds')

plot_grid(PWM_comp, stat_comp, signif_count, ncol = 3, 
          labels = c('A', 'B', 'C'), 
          label_fontface = 'plain', label_size = 14)-> stat_fig

ggsave(stat_fig, filename = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS7.svg', 
       width = 10, height = 5)


# supp fig validation screens
val_screen<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/validation_screen.rds')

spe_val<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/validation_specificity.rds')

max_val<-
  readRDS( '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/validation_max.rds')
val_avail<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/validation_avail.rds')

plot_grid(max_val+
            ylab('PPI signal\n(mean corrected AUC)')+
            theme(legend.position = 'none'), 
          spe_val+
            ylab('PPI signal\n(mean corrected AUC)'), val_screen, val_avail, 
          align = 'vh', axis = 'rltb')

plot_grid(max_val+
            ylab('PPI signal\n(mean corrected AUC)')+
            xlab('F[3]-peptide')+
            theme(legend.position = 'none', 
                  axis.text.x = element_text(vjust = 0.5)),
          val_screen+
            theme(legend.position = 'bottom'), nrow = 2, align = 'v', axis = 'l', 
          labels = c('A', 'C'), 
          label_fontface = 'plain', label_size = 14)->left


plot_grid( spe_val+
             ylab('PPI signal\n(mean corrected AUC)')+
             xlab('F[3]-peptide')+
             theme( axis.text.x = element_text(vjust = 0.5), 
                    legend.position = 'none'),
           val_avail, 
           align = 'v', axis = 'l', nrow =2,
           labels = c('B', 'D'), 
           label_fontface = 'plain', label_size = 14)->rigth


plot_grid(left, rigth, rel_widths = c(1, 0.75))-> SuppVal

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/Supp_validation_gc.svg',
       width = 14, height =8)


# Figure S : peptide properties

size_prop<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/size_prop.rds')
aa_prop<-readRDS( '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/aa_hydro.rds')+
  theme(legend.position = 'bottom')+
  guides(fill =guide_legend(title = 'peptide\nfamily'))

plot_grid(aa_prop, size_prop, 
          labels = c('A', 'B'), 
          label_fontface = 'plain', label_size = 14)->aa_prop_change

ggsave(aa_prop_change, filename = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/aa_properties_change.png', 
       height = 4, width = 6)
