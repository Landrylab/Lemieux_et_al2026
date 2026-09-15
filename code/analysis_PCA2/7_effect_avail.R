# AUTHOR : Pascale Lemieux
# Check various features of peptide sequences effect on availability

# import libraries and custom functions
library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# import median scores of the July2025 experiment (single mutant availability)
specificity <- 
  read_csv('signif_score_welch20_filter.csv')

avail_screen2 <- subset(specificity, pool == 'P1' & condition == 'MTX')

# import the scores from the February2025 experiment (double mutant availability)
avail <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/avail_signif_score.csv')


# select only single mutant from the first screen Feb2025
avail_screen1_single <- subset(avail, deg1 == max1 | deg2 == max2)

avail_screen1_single[avail_screen1_single$deg1 != avail_screen1_single$max1, 'pos'] <- 
  avail_screen1_single[avail_screen1_single$deg1 != avail_screen1_single$max1, 'pos1']

avail_screen1_single[avail_screen1_single$deg2 != avail_screen1_single$max2, 'pos'] <- 
  avail_screen1_single[avail_screen1_single$deg2 != avail_screen1_single$max2, 'pos2'] 

avail_screen1_single[avail_screen1_single$deg1 != avail_screen1_single$max1, 'deg'] <- 
  avail_screen1_single[avail_screen1_single$deg1 != avail_screen1_single$max1, 'deg1']

avail_screen1_single[avail_screen1_single$deg2 != avail_screen1_single$max2, 'deg'] <-
  avail_screen1_single[avail_screen1_single$deg2 != avail_screen1_single$max2, 'deg2']

avail_screen1_single%<>%
  subset(family != 'reference')%>%
  mutate(family = as.numeric(family))


avail_screen1_single$pos <- 
  as.numeric(avail_screen1_single$pos)

avail_screen1_single%<>%
    select(aa_seq, max_pep, family, med_norm_sel_coef, availability, pos, deg)

colnames(avail_screen1_single)[c(2,4,5)]<-c('peptide_max', 'med_norm', 'side')


# select unique aa sequences from the July 2025 screen
avail_screen2%>%
select(aa_seq, family, condition, med_norm, med_sel, pos, deg, peptide_max, side)%>%
  subset(family != 'reference')%>%
  dplyr::mutate(family= as.numeric(family))%>%
  unique()-> avail_screen2_u




# combine both experiment for single mutants
all_avail_single <- 
bind_rows(avail_screen2_u, 
          avail_screen1_single)

all_avail_single$condition <- 
  replace_na(all_avail_single$condition, 'MTX')

# order aa based on their chemical properties
aa <- c('R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'A', 'V', 'I','L', 'M', 'F', 'Y', 'W', 'C', 'G', 'P', '*' )

all_avail_single$deg <- 
factor(all_avail_single$deg, 
       levels= aa)

# Visualize availability for all the single mutants from both experiments

all_avail_single%>%
  dplyr::filter(family %in% c('363', '366', '385'))->sub_df


avail_position<-
ggplot(sub_df)+
  facet_grid(cols = vars(family),
             scales ='free')+
  geom_tile(aes(x = as.factor(pos), y =deg, fill = med_norm))+
  geom_point(data = subset(sub_df, peptide_max),
             aes(x=pos, y=deg), shape = 4, color ='red')+
  xlab('position')+
  ylab('amino acid')+
  scale_y_discrete(limits = aa, na.translate = F)+
  scale_x_discrete(na.translate = F)+
  scale_fill_viridis_c(option = 'E', na.value = 'white',  limits =c(0, 1), oob = scales::oob_squish, 
  name ='Avail. score')+
  t+
  theme(legend.position ='bottom')


ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/position_effect_avail.svg',avail_position, width = 3*2.5, height = 4)
saveRDS(avail_position, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/position_effect_avail.rds')


# compute feature (median, standard deviation) per position accross peptide for stop codons
all_avail_single_stop <- all_avail_single[all_avail_single$deg == '*', ] %>%
  dplyr::group_by(pos) %>%
  dplyr::summarise(
    med_pos = median(med_norm, na.rm = TRUE),
    sd_pos = sd(med_norm, na.rm = TRUE),
    n = n(),
    se_pos = sd_pos / sqrt(n)
  )

# plot the trend of position on availability for stop codons

avail_stop<-
ggplot(all_avail_single_stop, aes(x = as.factor(pos), y = med_pos)) +
  geom_smooth(aes(x = as.numeric(as.character(pos)), y = med_pos), method = 'lm', color = 'black', linetype = 'dashed', inherit.aes = FALSE)+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = med_pos - se_pos, ymax = med_pos + se_pos), width = 0.2) +
  scale_x_discrete(na.translate = FALSE) +
  stat_cor(aes(x = as.numeric(as.character(pos)), y = med_pos), 
  method = 'spearman', label.x = 5, label.y = 0.9, size = 3, inherit.aes = FALSE, cor.coef.name = 'rho', 
  label.sep = '\n', ) +
  xlab('stop codon position') +
  ylab('Avail. score') +
  t

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/med_stop_position_effect_avail.png', width = 6, height = 4)
saveRDS(avail_stop, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/med_stop_position_effect_avail.rds')

# select sequences for stop codon validation effect on availability
all_avail_single%>%
  filter(deg == '*' & pos %in% c(2,5,8,10) & family %in% c(363, 366, 385))->valid_stop
  
valid_stop%>%
  dplyr::summarise(max = max(med_norm), 
            min = min(med_norm)) 

# Compare availability and peptide stickiness (based on Levy 2014 and update paper 2022)

# import stickiness score for each amino acid
sticky <- read_table('~/PL_projects/PL_papers/PPI_optimization_paper/data/sticky_score.txt')
sticky2022<- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/stickiness2022.csv')
sticky_IDP<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/stickiness_IDP.csv')[, -1]

hydrophobicity<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/hydrophobicity.csv')
colnames(hydrophobicity)[1]<-'AA'

sticky%<>%
  dplyr::arrange(sticky_score)%>%
  dplyr::mutate(AA = factor(AA, levels = AA, ordered = T))


hydrophobicity%<>%
  dplyr::arrange(Cornette)%>%
  dplyr::mutate(AA = factor(AA, levels = AA, ordered = T))

my_gradient <- colorRampPalette(c('#373737', '#fde333ff'))

ggplot(hydrophobicity)+
  geom_col(aes(x=AA, y=Cornette, fill = AA))+
  scale_fill_manual(values = my_gradient(20))+
  #scale_x_discrete(limits = AA)+
  t+
  ylab('Cornette hydrophobicity')+
  theme(legend.position = 'none')->stiky_scale

saveRDS(stiky_scale, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/hydro_scale.rds')


colnames(sticky_IDP)[c(1,3)]<-c('AA', 'sticky_score')

comp_sticky<-
  full_join(sticky, sticky2022, join_by('AA'), 
            suffix = c('2012', '2022'))

comp_sticky<-
  full_join(comp_sticky, sticky_IDP[, c(1,3)], join_by(AA), 
            suffix = c('', '.IDP'))

colnames(comp_sticky)[4]<-'sticky_scoreIDP'

comp_sticky<-
full_join(comp_sticky, hydrophobicity, join_by(AA))


ggplot(comp_sticky, aes(x =sticky_score2022, y = Cornette))+
  geom_point()+
  stat_cor()+
  t

Cornette<-
  comp_sticky%>%
  select('AA', 'Cornette')


colnames(Cornette)[2]<-'sticky_score'

# select all aa sequence from both screens
avail%>%
  select(aa_seq, family, med_norm_sel_coef)%>%
  unique() ->avail_screen1


colnames(avail_screen1)[3]<-'med_norm'

avail_screen2%>%
  select(aa_seq, family, med_norm)%>%
  unique() -> avail_screen2

# combine all sequences 
all_Avail <- 
  bind_rows(avail_screen1, avail_screen2, .id ='screen')

# identify stop codons
stop <- 
  grepl(all_Avail$aa_seq, pattern = '*', fixed = T)

# modify the aa sequence based on the stop codon
aa_seq <- 
  strsplit(all_Avail$aa_seq, '*', fixed = T)
aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')

# compute the sticky score for each aa sequence 

for(j in c('sticky', 'sticky2022', 'sticky_IDP', 'Cornette')){

  for (i in 1:length(aa_seq)) {
  if (length(aa_seq[[i]]) == 0){
    all_Avail[i, j] <- 0
    next
  }
  all_Avail[i, j] <- get_stickyness(aa_seq[[i]], get(j))
  }
}  

all_Avail$family <- 
  as.factor(all_Avail$family)
all_Avail$stop <- stop

# Plot sticky score vs availability separated for each screen 
# 1: double mutants from Feb2025, 2: single mutants from July2025

all_Avail%>%
  filter(family != 'reference' )%>%
ggplot()+
  facet_grid(cols = vars(screen))+
  geom_point(aes(x=med_norm, y = sticky, color = family, shape = stop))+
  scale_shape_manual(values = c(1,4))+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey30'))+
  stat_cor(aes(x=med_norm, y = sticky), method = 'spearman', cor.coef.name = 'rho', label.y = -2, label.sep = '\n', size=3)+
  xlab('availability score')+
  ylab('stickiness score\n(Levy et al. 2012)')+
  t+
  theme(legend.position = 'bottom')->sticky2012


all_Avail%>%
  filter(family != 'reference' )%>%
  ggplot()+
  facet_grid(cols = vars(screen))+
  geom_point(aes(x=med_norm, y = sticky2022, color = family, shape = stop))+
  scale_shape_manual(values = c(1,4))+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey30'))+
  stat_cor(aes(x=med_norm, y = sticky2022), method = 'spearman', cor.coef.name = 'rho', label.y = 0, label.sep = '\n',  size=3)+
  xlab('availability score')+
  ylab('stickiness score\n(Villegas and Levy, 2022)')+
  t+
  theme(legend.position = 'bottom')->Supp_sticky_2022


all_Avail%>%
  filter(family != 'reference' & !stop)%>%
  ggplot()+
  facet_grid(cols = vars(screen))+
  geom_point(aes(x=med_norm, y = sticky_IDP, color = family, shape = stop))+
  scale_shape_manual(values = c(1,4))+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey30'))+
  stat_cor(aes(x=med_norm, y = sticky_IDP), method = 'spearman', cor.coef.name = 'rho', label.y = 2, label.sep = '\n',  size=3)+
  xlab('availability score')+
  ylab('stickiness score\n(Cao et al., 2026)')+
  t+
  theme(legend.position = 'bottom')+
  guides(color = guide_legend(title = 'peptide family'), 
         shape = guide_legend(nrow=2))->Supp_stiky_IDP


all_Avail%>%
  filter(family != 'reference'& !stop & screen==1)%>%
  ggplot()+
  facet_grid(
             cols = vars(family))+
  geom_point(aes(x=Cornette, y = med_norm, color = family), shape=1)+
  #scale_shape_manual(values = c(1,4))+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey30'))+
  stat_cor(aes(x=med_norm, y = Cornette), method = 'spearman', cor.coef.name = 'rho', label.y = 0.9, label.x = 10, label.sep = '\n',  size=3)+
  ylab('availability score')+
  xlab('Hydrophobicity (Cornell scale)')+
  t+
  theme(legend.position = 'none')->Supp_hydro

#+
  guides(color = guide_legend(title = 'peptide family'), 
         shape = guide_legend(nrow=2))->Supp_hydro


  all_Avail%>%
    filter(family != 'reference'& screen==1)%>%
    ggplot()+
    geom_point(aes(x=Cornette, y = med_norm, color = family, shape = stop))+
    scale_shape_manual(values = c(1,4))+
    scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF", 'grey30'))+
    stat_cor(aes(x=med_norm, y = Cornette), method = 'spearman', cor.coef.name = 'rho', label.y = 0.9, label.x = 10, label.sep = '\n',  size=3)+
    ylab('availability score')+
    xlab('Hydrophobicity (Cornell scale)')+
    t+
    theme(legend.position = 'none')->Supp_hydro_c
  

  leg_stik<-get_legend(Supp_stiky_IDP)


ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/avail_vs_stickiness.png', stiky2012, width = 10, height = 6)
saveRDS(sticky2012, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/avail_vs_stickiness.rds')
saveRDS(Supp_hydro, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/avail_vs_hydro_all.rds')
saveRDS(Supp_hydro_c, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/avail_vs_hydro.rds')


plot_grid(Supp_sticky_2022+theme(legend.position = 'none'),
          Supp_stiky_IDP+theme(legend.position = 'none'), ncol =2, 
          labels = c('A', 'B'), label_size = 14, label_fontface = 'plain')->sticki_plot

plot_grid(sticki_plot, leg_stik, rel_heights = c(1, 0.2), nrow = 2)->SuppSticki

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/SuppStiki.png', SuppSticki, 
       width=8, height = 5)


# save the validation peptide with stop codon and stickiness score
left_join(valid_stop, all_Avail[, c('aa_seq', 'screen', 'sticky_score')])-> valid_stop
real_seq<-strsplit(valid_stop$aa_seq, split ='*', fixed = TRUE)

real_seq<-
  unlist(lapply(real_seq, `[[`, 1))

valid_stop$express_seq<-real_seq
write_csv(valid_stop, 'peptide_validation_stop.csv')

# compare availability vs fitness 
# miss the dna sequence for the merging between fitness and avail.s

# import fitness selection coefficient

# avail score from July2025
avail_fitness_screen2<-
  specificity%>%
  dplyr::filter(condition=='DMSO' & pool == 'P2')

# import the avail score from Feb2025

avail_fitness_screen2%>%
  #select(aa_seq, sticky_score)%>%
  unique()%>%
  left_join(all_Avail, join_by('aa_seq', 'family'))->avail_vs_fitness


# plot fitness vs availability
avail_vs_fitness%>%
dplyr::filter(family != 'reference')%>%
ggplot()+
geom_point(aes(x = med_norm.y, y = med_sel, color= family, shape = stop))+
stat_cor(aes(x = med_norm.y, y = med_sel), method = 'spearman', cor.coef.name = 'rho')+
scale_shape_manual(values = c(1,4))+
scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
labs(x = 'Avail. score', y = 'selection coefficient')+
t+
theme(legend.position = 'bottom')+
guides(shape = 'none')->fit_vs_avail
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/fit_vs_avail.png', fit_vs_avail, height = 5, width = 5)



## Recompute stickiness for the missing sequences
# identify stop codons
stop <- 
  grepl(avail_vs_fitness$aa_seq, pattern = '*', fixed = T)

# modify the aa sequence based on the stop codon
aa_seq <- 
  strsplit(avail_vs_fitness$aa_seq, '*', fixed = T)
aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')

# compute the stiky score for each aa sequence 
for (i in 1:length(aa_seq)) {
  if (length(aa_seq[[i]]) == 0){
    avail_vs_fitness[i, 'sticky_score'] <- 0
    next
  }
  avail_vs_fitness[i, 'sticky_score'] <- get_stickyness(aa_seq[[i]], sticky)
}


avail_vs_fitness%>%
  subset(family != 'reference')%>%
  ggplot()+
  geom_point(aes(x = sticky_score, y = med_sel, color= family, shape = stop))+
  stat_cor(aes(x = sticky_score, y = med_sel), method = 'spearman', cor.coef.name = 'rho')+
  scale_shape_manual(values = c(1,4))+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  labs(x = 'stickiness score', y = 'selection coefficient')+
  t+
  theme(legend.position = 'bottom')+
  guides(shape = 'none')->fit_vs_sticky

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/fit_vs_sticky.png', fit_vs_sticky, height = 5, width = 5)



