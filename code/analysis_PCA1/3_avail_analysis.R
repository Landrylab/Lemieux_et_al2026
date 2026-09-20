# AUTHOR : Pascale Lemieux
### Availability assay analysis ###
#import packages and functions

library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggsci)


source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/data/PCA1/')

# import selection coefficient
long_delta_avail <- read_csv('sel_coef_availability.csv')

# scale avail signal at T2 using sh3.1 and median between 304_1 & 250_4

long_delta_avail%>%
  filter(family %in% c('250_4', 'pdz.2') & timepoint_end == 'T2')%>%
  select(replicate, family, sel_coef) -> ref_avail

pivot_wider(ref_avail, values_from = sel_coef, names_from = family)-> ref_avail

colnames(ref_avail)[2:3] <- c('neg_avail', 'pos_avail')

# merge with ref values
long_delta_avail_T2 <- 
  merge(long_delta_avail[long_delta_avail$timepoint_end == 'T2', ], 
        ref_avail, 
        by = 'replicate')

# compute scaling to ref values
long_delta_avail_T2 %<>%
  mutate(norm_sel_coef = (sel_coef - neg_avail)/(pos_avail - neg_avail))


before_norm_T2 <- 
  ggplot(long_delta_avail_T2)+
  #facet_wrap(vars(timepoint_end))+
  geom_density(aes(sel_coef, color =replicate))+
  scale_color_cosmic()+
  xlab('selection coefficient (s)')+
  t

post_norm_T2 <- 
  ggplot(long_delta_avail_T2)+
  #facet_wrap(vars(timepoint_end))+
  geom_density(aes(norm_sel_coef, color =replicate))+
  xlab('normalized s')+
  scale_color_cosmic()+t+theme(legend.position = 'bottom')


legend = get_plot_component(post_norm_T2, 'guide-box-bottom', return_all = TRUE)
saveRDS(legend, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/legend_FigS4BC.rds')


FigS4B <- 
plot_grid(before_norm_T2+theme(legend.position = 'none'), post_norm_T2+theme(legend.position = 'none'), nrow = 2, 
          align = 'v', axis = 'lr', rel_heights = c(1, 1))



saveRDS(FigS4B, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4B.rds')

# check if reference becomes fitter with time 
# scale T1 also with reference sequences

long_delta_avail%>%
  filter(family %in% c('250_4', 'pdz.2') & timepoint_end == 'T1')%>%
  select(replicate, family, sel_coef) -> ref_avail

pivot_wider(ref_avail, values_from = sel_coef, names_from = family)-> ref_avail

colnames(ref_avail)[2:3] <- c('neg_avail',  'pos_avail')


# merge with ref values
long_delta_avail_T1 <- 
  merge(long_delta_avail[long_delta_avail$timepoint_end == 'T1', ], 
        ref_avail, 
        by = 'replicate')

# compute scaling to ref values
long_delta_avail_T1 %<>%
  mutate(norm_sel_coef = (sel_coef - neg_avail)/(pos_avail - neg_avail))

comp_timepoint <- bind_rows(long_delta_avail_T1, long_delta_avail_T2)

norm_coef_rep <- 
ggplot(comp_timepoint)+
  geom_violin(aes(y = norm_sel_coef, x = replicate, fill =timepoint_end))+
  scale_fill_grey()+
  ylab('normalized s')+
  t

coef_rep <- 
  ggplot(comp_timepoint)+
  ylab('selection coefficient (s)')+
  geom_violin(aes(y = sel_coef, x = replicate, fill =timepoint_end))+
  scale_fill_grey()+
  theme_classic2()+
  t

FigS4C <- 
plot_grid(coef_rep+theme(legend.position = 'none'), norm_coef_rep+theme(legend.position = 'bottom'), nrow = 2, 
          align = 'v', axis = 'lr', rel_heights = c(1, 1.2))


saveRDS(FigS4C, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4C.rds')

## norm sel coefficient is constant across timepoints and triplicates

long_delta_avail_T2$max_pep <- 
  c(long_delta_avail_T2$deg1 ==long_delta_avail_T2$max1 &long_delta_avail_T2$deg2 ==long_delta_avail_T2$max2)


## Check replicates 
replicate_check_avail<- 
  pivot_wider(long_delta_avail_T2[, c(1:13,22,23)], 
              names_from = replicate, 
              values_from = norm_sel_coef)

replicate_check_avail[!grepl(pattern = 'PRM_', replicate_check_avail$family), 'ref_family'] <- TRUE
replicate_check_avail$ref_family <- replace_na(replicate_check_avail$ref_family, replace = FALSE)

replicate_check_avail[replicate_check_avail$ref_family, 'alpha_val'] <- 1
replicate_check_avail[!replicate_check_avail$ref_family, 'alpha_val'] <- 0.1

# Use custom function to generate correlation graph between replicates
Fig3A <- 
  ggpairs(replicate_check_avail, columns = 14:16,
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

library(gridExtra)
library(gridGraphics)
library(patchwork)
library(ggplotify)
saveRDS(Fig3A, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig3A.rds')
Fig3A<- grid.grabExpr(grid.draw(Fig3A))



# Compare with growth curves measurements

gc_data <- 
read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/condensed_gc_PPI_Avail.csv')[,-1]

ref_sequences <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ref_sequences.csv')

gc_data%<>%
  filter(strain =='PL17')

ref_sequences <- 
merge(y = gc_data, 
      x = ref_sequences, 
      by.y = 'D3', 
      by.x = 'name')

ref_comp <- 
  tibble(
  merge(y = long_delta_avail_T2, 
        x = ref_sequences, 
        by.y = c('aa_seq', 'sequence'), 
        by.x = c('peptide', 'trim')))

# simplify ref_comp df
ref_comp <- ref_comp[, c(1:3, 8,9, 10,12:21,26, 29,30)]


# change between 'sequence' and 'aa_seq' for dna vs peptide unique sequence 
# not better with dna sequences
sel_coef_summary <- 
  data_summary(long_delta_avail_T2, 'sel_coef', c('aa_seq'))

norm_sel_coef_summary <- 
  data_summary(long_delta_avail_T2, 'norm_sel_coef', c('aa_seq'))

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

colnames(sd_gc)[1] <- 'sequence'
comp_w_gc <- 
inner_join(sd_gc, comp_w_gc, by = "sequence")


# auc vs norm_sel_coef 
# remove the sequence with a strong outlier
Fig3B <- 
ggplot(data = comp_w_gc[comp_w_gc$sequence != 'GLLWSRLAPSVPVL', ],
       aes(x = auc,y = norm_sel_coef)) + 
  geom_pointrange(aes(xmin = auc-sd_auc,xmax = auc+sd_auc), color = 'darkgrey', size = 0.1, linewidth = 0.8) + 
  geom_pointrange(aes(ymin = norm_sel_coef-sd_norm_sel_coef, ymax = norm_sel_coef+sd_norm_sel_coef), 
                  color = 'darkgrey', size = 0.1, linewidth = 0.8)+
  geom_point(size =1, color = 'black') + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 3, label.sep = '\n', 
           label.x = 5, label.y = 0.5)+
  ylab('Avail. score')+
  xlab('corrected AUC')+
  t

saveRDS(Fig3B, file = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig3B.rds')

#  Generate heatmap

long_delta_avail_T2[!grepl(pattern = 'PRM', long_delta_avail_T2$family), 'family'] <- 'reference'
long_delta_avail_T2[long_delta_avail_T2$family == 'reference', c('deg1', 'deg2', 'pos1', 'pos2', 'max1', 'max2')] <- 
  long_delta_avail_T2[long_delta_avail_T2$family == 'reference', 'codon1']

# compute median per aa sequence

long_delta_avail_T2[, c(4,7:13, 22)]%>%
  dplyr::group_by(aa_seq)%>%
  dplyr::mutate(med_norm = median(norm_sel_coef))%>%
  select(aa_seq, deg1, deg2, pos1, pos2, max1, max2, family, med_norm)%>%
  unique()-> med_norm


max_pep_sel <- 
  med_norm[med_norm$deg1 == med_norm$max1 & med_norm$deg2 == med_norm$max2 , ]

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
gsub(med_norm$family, pattern = 'PRM_0', replacement = '')


max_pep_sel$name <-
  gsub(max_pep_sel$family, pattern = 'PRM_0', replacement = '')

FigS6 <- 
ggplot(med_norm[med_norm$name != 'reference', ])+
  facet_wrap(vars(name), scales = 'free')+
  geom_tile(aes(deg1_order, deg2_order, fill = med_norm))+
  geom_point(data = max_pep_sel[max_pep_sel$name != 'reference', ], aes(deg1, deg2), 
             color = 'red', shape = 4)+
  scale_fill_viridis_c(option = 'E')+
  t+
  xlab('deg. position 1')+
  ylab('deg. position 2')+
  theme(legend.position.inside = c(0.6, 0.25))+
  guides(fill = guide_colorbar(title = 'availability\nscore', 
                               position = 'inside',
                               direction = 'horizontal',
                               
                               theme = theme(
                                 legend.key.width  = unit(10, "lines"),
                                 legend.key.height = unit(1.5, "lines")
                               )))

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS15.png', FigS6, width = 8, height = 8)

# Distribution of effect
Fig3C <- 
ggplot(med_norm[med_norm$name !=  'reference', ])+
  geom_violin(aes(med_norm, name), trim = TRUE, width = 0.8, draw_quantiles = c(0.25,0.75), color = 'black', fill = 'transparent')+
  geom_point(data = max_pep_sel[(max_pep_sel$name !=  'reference'), ], 
             aes(med_norm, name, color = name), size=2)+
  geom_vline(aes(xintercept=0), linetype = 'dashed')+
  stat_summary(aes(x = med_norm, name), fun.y=median, geom="point", size=2, color = 'black', shape = 4)+
  scale_color_viridis_d(option = 'G')+
  geom_vline(aes(xintercept=1), linetype = 'dashed')+
  xlab('Avail. score')+
  ylab('peptide family')+
  t+
  theme(legend.position = 'none')

saveRDS(Fig3C, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigXC.rds')

# save median normalized score for further analysis

write_csv(long_delta_avail_T2, 'norm_avail_sel_coef.csv')
write_csv(med_norm, 'med_avail_score.csv')

