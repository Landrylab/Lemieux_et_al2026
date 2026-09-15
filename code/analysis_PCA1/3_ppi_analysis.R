#AUTHOR: Pascale Lemieux
### PPI assay analysis ###
#import packages

library(tidyverse)
library(Biostrings)
library(magrittr)
library(ggpubr)
library(cowplot)
library(ggsci)
library(GGally)
library(grid)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/')

# import selection coefficient
long_delta_ppi <- read_csv('./data/PCA1/sel_coef_ppi.csv')

# scale PPI signal at T2 using sh3.1 and 250_3

long_delta_ppi%>%
  dplyr::filter(family %in% c("250_3", "sh3.1") & timepoint_end == 'T2')%>%
  dplyr::select(replicate, family, sel_coef) -> ref_PPI

pivot_wider(ref_PPI, values_from = sel_coef, names_from = family)-> ref_PPI

colnames(ref_PPI)[2:3] <- c('neg_PPI', 'pos_PPI')

# merge with ref values
long_delta_ppi_T2 <- 
  merge(long_delta_ppi[long_delta_ppi$timepoint_end == 'T2', ], 
        ref_PPI, 
        by = 'replicate')

# compute scaling to ref values
long_delta_ppi_T2 %<>%
  mutate(norm_sel_coef = (sel_coef - neg_PPI)/(pos_PPI - neg_PPI))

before_norm_T2 <- 
ggplot(long_delta_ppi_T2)+
  #facet_wrap(vars(timepoint_end))+
  geom_density(aes(sel_coef, color =replicate))+
  scale_color_cosmic()+
  xlab('selection coefficient (s)')+
  t

post_norm_T2 <- 
  ggplot(long_delta_ppi_T2)+
  #facet_wrap(vars(timepoint_end))+
  geom_density(aes(norm_sel_coef, color =replicate))+
  scale_color_cosmic()+
  xlab('normalized s')+
  t+
  theme(legend.position = 'bottom')

FigS4D <- 
plot_grid(before_norm_T2+theme(legend.position = 'none'), post_norm_T2+theme(legend.position = 'none'), nrow = 2, 
          rel_heights = c(1, 1))

saveRDS(FigS4D, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4E.rds')
# check if reference becomes fitter with time 
# scale T1 also with reference sequences

long_delta_ppi%>%
  dplyr::filter(family %in% c("250_3", "sh3.1") & timepoint_end == 'T1')%>%
  dplyr::select(replicate, family, sel_coef) -> ref_PPI

pivot_wider(ref_PPI, values_from = sel_coef, names_from = family)-> ref_PPI

colnames(ref_PPI)[2:3] <- c('neg_PPI','pos_PPI')


# merge with ref values
long_delta_ppi_T1 <- 
  merge(long_delta_ppi[long_delta_ppi$timepoint_end == 'T1', ], 
        ref_PPI, 
        by = 'replicate')

# compute scaling to ref values
long_delta_ppi_T1 %<>%
  mutate(norm_sel_coef = (sel_coef - neg_PPI)/(pos_PPI - neg_PPI))


comp_timepoint <- bind_rows(long_delta_ppi_T1, long_delta_ppi_T2)

FigS4E <- 
ggplot(comp_timepoint)+
  geom_violin(aes(y = norm_sel_coef, x = replicate, fill =timepoint_end))+
  ylab('normalized s')+
  scale_fill_grey()+
  t+
  theme(legend.position = 'bottom')+
  guides(fill = guide_legend(title = 'end timepoint'))

saveRDS(FigS4E, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4F.rds')

## norm sel coefficient is constant across timepoints and triplicates
long_delta_ppi_T2$max_pep <- 
  c(long_delta_ppi_T2$deg1 ==long_delta_ppi_T2$max1 &long_delta_ppi_T2$deg2 ==long_delta_ppi_T2$max2)


## Check replicates post normalization
replicate_check_PPI<- 
  pivot_wider(long_delta_ppi_T2[, c(1:13, 22,23)], 
              names_from = replicate, 
              values_from = norm_sel_coef)

replicate_check_PPI[!grepl(pattern = 'PRM_', replicate_check_PPI$family), 'ref_family'] <- TRUE
replicate_check_PPI$ref_family <- replace_na(replicate_check_PPI$ref_family, replace = FALSE)

replicate_check_PPI[replicate_check_PPI$ref_family, 'alpha_val'] <- 1
replicate_check_PPI[!replicate_check_PPI$ref_family, 'alpha_val'] <- 0.1


Fig2A <- 
ggpairs(replicate_check_PPI, columns = 14:16,
        # Replace with your actual replicate columns
        mapping = aes(color = ref_family, alpha = alpha_val),
        upper = list(continuous = my_cor),
        lower = list(continuous = my_scatter_double_layer),
        diag  = list(continuous = my_diag))+
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 20), unit = 'mm') # Left margin)


Fig2A<- grid.grabExpr(grid.draw(Fig2A))

saveRDS(Fig2A, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig2A.rds')

## save normalized sel_coef
write_csv(long_delta_ppi_T2, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/norm_sel_coef_ppi.csv')

long_delta_ppi_T2 <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/norm_sel_coef_ppi.csv')
#top peptides are also toward the higher end of PPI strength

# Compare with growth curves measurements

gc_data <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/condensed_gc_PPI_Avail.csv')[,-1]

gc_data%<>%
  dplyr::filter(strain =='BY4741' & D12 != 'empty')

# remove bad gc data
d3 <- gsub(pattern = 'PRM_0', replacement = '', gc_data$D3)

d3 <- 
  gsub(pattern = '_.', replacement = '', d3)

gc_data <- 
  gc_data[gc_data$D12 == d3 | gc_data$D12 %in% c('PDZ', 'SH3'), ]

gc_data <- 
gc_data[!(gc_data$D3 == 'PRM_0363_4' & gc_data$plate == 7), ]

gc_data <- 
  gc_data[!(gc_data$D3 == 'PRM_0299_5' & gc_data$plate == 5), ]

# no sh31 and SH3 in reference sequence , to add 
ref_sequences <- read_csv('~/PL_projects/sequencing/Screen_PPI_Avail_Feb2025/ref_sequences.csv')

ref_sequences <- 
bind_rows(ref_sequences,
tibble(name = 'sh31', trim = 'CCTCCTCCTGCCCTACCCCCAAAAAGGAGACGTTAA', peptide = 'PPPALPPKRRR'))

ref_sequences <- 
  merge(y = gc_data, 
        x = ref_sequences, 
        by.y = 'D3', 
        by.x = 'name')


# change between 'sequence' and 'aa_seq' for dna vs peptide unique sequence 
# not better with dna sequences
sel_coef_summary <- 
  data_summary(long_delta_ppi_T2, 'sel_coef', c('aa_seq'))

norm_sel_coef_summary <- 
  data_summary(long_delta_ppi_T2, 'norm_sel_coef', c('aa_seq'))

auc_summary <- 
  data_summary(ref_sequences, 'auc', c('peptide'))

dgr_summary <- 
  data_summary(ref_sequences, 'dgr', c('peptide'))


colnames(sel_coef_summary)[3] <- 'sd_sel_coef'
colnames(norm_sel_coef_summary)[3] <- 'sd_norm_sel_coef'

colnames(auc_summary)[3] <- 'sd_auc'
colnames(dgr_summary)[3] <- 'sd_dgr'

comp_w_gc <- 
  full_join(sel_coef_summary, norm_sel_coef_summary)

colnames(comp_w_gc)[1] <- 'sequence'

sd_gc <- 
  full_join(auc_summary, dgr_summary)

# comparison between dgr and auc -> auc has less background than dgr
ggplot(data = sd_gc,aes(x = auc,y = dgr)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = auc-sd_auc,xmax = auc+sd_auc)) + 
  geom_errorbar(aes(ymin = dgr-sd_dgr,ymax = dgr+sd_dgr))+
  theme_classic2()

colnames(sd_gc)[1] <- 'sequence'

comp_w_gc <- 
  inner_join(sd_gc, comp_w_gc, by = 'sequence')

# auc vs norm_sel_coef
# remove sd outlier WFPLVDYLVHIS

Fig2B <- 
ggplot(data = comp_w_gc[comp_w_gc$sequence != 'WFPLVDYLVHIS', ],
       aes(x = auc,y = norm_sel_coef)) + 
  geom_pointrange(aes(xmin = auc-sd_auc,xmax = auc+sd_auc), color = 'darkgrey', size = 0.1, linewidth = 0.8) + 
  geom_pointrange(aes(ymin = norm_sel_coef-sd_norm_sel_coef, ymax = norm_sel_coef+sd_norm_sel_coef), 
                color = 'darkgrey', size = 0.1, linewidth = 0.8)+
  geom_point(size =1, color = 'black') + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 3, label.sep = '\n', label.x = 10, label.y = 0.8)+
  ylab('PPI score')+
  xlab('corrected AUC')+
  t


library(gridExtra)
library(gridGraphics)
library(patchwork)
library(ggplotify)

#Fig2C<- grid.grabExpr(grid.draw(Fig2A))
Fig2B <- ggplotGrob(Fig2B)

saveRDS(Fig2A, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig2A.rds')
saveRDS(Fig2B, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig2B.rds')

# Generate heatmap for aa 

long_delta_ppi_T2[!grepl(pattern = 'PRM', long_delta_ppi_T2$family), 'family'] <- 'reference'
long_delta_ppi_T2[long_delta_ppi_T2$family == 'reference', c('deg1', 'deg2', 'pos1', 'pos2', 'max1', 'max2')] <- 
  long_delta_ppi_T2[long_delta_ppi_T2$family == 'reference', 'codon1']

# compute median per aa sequence

long_delta_ppi_T2[, c(4,7:13, 22)]%>%
  dplyr::group_by(aa_seq)%>%
  dplyr::mutate(med_norm = median(norm_sel_coef))%>%
  dplyr::select(aa_seq, deg1, deg2, pos1, pos2, max1, max2, family, med_norm)%>%
  unique()-> med_norm


max_pep_sel <- 
  med_norm[med_norm$deg1 == med_norm$max1 & med_norm$deg2 == med_norm$max2 , ]%>%
  drop_na()

#order by amino acid per chemical properties
# positive : R, H, K
# negative : D, E
# Polar : S, T, N, Q 
# hydrophobic small: A, V, I, L, 
# hydrophobic large : M, F, Y, W
# special : C, G, P

aa <- 
c('R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'A', 'V', 'I','L', 'M', 'F', 'Y', 'W', 'C', 'G', 'P', '*' )

med_norm$deg1_order <- 
factor(med_norm$deg1, 
       labels = aa, 
       levels = aa)

med_norm$deg2_order <- 
  factor(med_norm$deg2, 
         labels = aa, 
         levels = aa)

med_norm$name <- 
gsub(med_norm$family, pattern = 'PRM_0', replacement= '')
max_pep_sel$name <- 
gsub(max_pep_sel$family, pattern = 'PRM_0', replacement= '')

#FigS4A <- 
med_norm[!(med_norm$name %in% c('reference', '366')), ]%>%
ggplot()+
  facet_wrap(vars(name), scales = 'free')+
  geom_tile(aes(deg1_order, deg2_order, fill = med_norm), linejoin = 'bevel')+
  geom_point(data = max_pep_sel[!(max_pep_sel$name %in% c('reference', '366')), ], 
             aes(deg1, deg2), color = 'black', shape = 4, size = 2.5)+
  scale_fill_viridis_c(option = 'F', direction = 1, end = 0.95)+
  xlab('deg. position 1')+
  ylab('deg. position 2')+
  t+
  theme(strip.text = element_text(color = 'black'))+
  guides(fill = guide_colorbar(title = 'PPI score', 
                               theme = theme(
                                 legend.key.width  = unit(1, "lines"),
                                 legend.key.height = unit(8, "lines")
                               )))

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS5.png', width = 9, height = 6)

med_norm%>%
  dplyr::filter(name == '366')%>%
  dplyr::group_by(deg1_order)%>%
  dplyr::summarise(mean_per_aa = median(med_norm))%>%
  dplyr::arrange(mean_per_aa)->aa_order1

aa_order1$deg1_order<-
factor(aa_order1$deg1_order, 
       levels = aa_order1$deg1_order)

med_norm$deg1_order<-
factor(med_norm$deg1_order, 
       levels = aa_order1$deg1_order)


med_norm%>%
  dplyr::filter(name == '366')%>%
  dplyr::group_by(deg2_order)%>%
  dplyr::summarise(mean_per_aa = median(med_norm))%>%
  dplyr::arrange(mean_per_aa)->aa_order2


aa_order2$deg2_order<-
  factor(aa_order2$deg2_order, 
         levels = aa_order2$deg2_order)


med_norm$deg2_order<-
  factor(med_norm$deg2_order, 
         levels = aa_order2$deg2_order)



Fig2E <-   
ggplot(med_norm[(med_norm$name == '366'), ])+
  facet_wrap(vars(name), scales = 'free')+
  geom_tile(aes(deg1_order, deg2_order, fill = med_norm), linejoin = 'bevel')+
  geom_point(data = max_pep_sel[(max_pep_sel$name ==  '366'), ], 
             aes(deg1, deg2), color = 'black', shape = 4, size = 2.5)+
  scale_fill_viridis_c(option = 'F', direction = 1, end = 0.95)+
  xlab('deg. position 1')+
  ylab('deg. position 2')+
  t+
  theme(strip.text = element_text(color = 'black'),
        legend.position = 'bottom')+
  guides(fill = guide_colorbar(title = 'PPI score', 
                               theme = theme(
                                 legend.key.width  = unit(8, "lines"),
                                 legend.key.height = unit(1, "lines")
                               )))
saveRDS(Fig2E, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig2E.rds')

x<-readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2E.rds')

med_norm$max <- 
  med_norm$aa_seq %in% max_pep_sel$aa_seq

med_norm <- 
left_join(med_norm, max_pep_sel[, c('name', 'med_norm')], by = c('name'), 
          suffix = c('', '.max'))

med_norm%<>%
  group_by(name)%>%
  mutate(higher = med_norm>=med_norm.max)


## after 4_ppi_stat_analysis.R
stat_signif<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ppi_signif_score.csv')


Fig2C <- 
ggplot(med_norm[med_norm$name !=  'reference', ])+
  geom_vline(aes(xintercept=0), linetype = 'dashed', color = 'grey45')+
   geom_jitter(data = stat_signif[stat_signif$binding=='stronger', ], 
             aes(med_norm_sel_coef, family), size=1.5, height = 0.2, color = 'grey45', shape =1, alpha=0.6)+
  geom_violin(aes(med_norm, name), trim = TRUE, width = 0.8, draw_quantiles = c(0.25,0.75), color = 'black', fill = 'transparent')+
  geom_point(data = max_pep_sel[(max_pep_sel$name !=  'reference'), ], 
             aes(med_norm, name, color = name), size=2)+
  stat_summary(aes(x = med_norm, name), fun=median, geom="point", size=2, color = 'black', shape = 4)+
  scale_color_viridis_d(option = 'G')+
  geom_vline(aes(xintercept=1), linetype = 'dashed', color = 'grey45')+
  xlab('PPI score')+
  ylab('peptide family')+
  t+
  theme(legend.position = 'none')


saveRDS(Fig2C, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2C.rds')


# save med norm score for further analysis
write_csv(med_norm, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/med_ppi_score.csv')






