# AUTHOR : Pascale Lemieux
### Selection coefficient analysis ###
#import packages and functions

library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggnewscale)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# import selection coefficient
norm_sel_coef <- read_csv('filtered_coefficient.csv')%>%
  select(condition, aa_seq, family, timepoint_end, replicate,sequence, pool, assay, sel_coef, pos1, pos2, pos, max1, max2, max, deg1, deg2, deg, double_mutant, norm_sel_coef)

table(norm_sel_coef$replicate)
table(norm_sel_coef$pool)

# Identify the reference peptide (maximal binder according to phage display)
norm_sel_coef%>%
subset((deg1 == max1 & deg2 == max2 )| 
        (deg == max ))%>%
        select(aa_seq)%>%
        unique()-> max_pep

norm_sel_coef[norm_sel_coef$aa_seq %in% max_pep$aa_seq, 'peptide_max']<-TRUE

# duplicate the values for the reference peptide to add a reference for all the position of the peptide family
norm_sel_coef%>%
  ungroup()%>%
  subset(peptide_max)%>%
  dplyr::mutate(pos = list(1:10)) %>%
  unnest(pos)%>%
  dplyr::mutate(
    max = map2_chr(aa_seq, pos, ~ substr(.x, .y, .y)),
    deg = map2_chr(aa_seq, pos, ~ substr(.x, .y, .y))
  )->x

norm_sel_coef%<>%
  bind_rows(x) 

# save the sel coefficient with the added reference peptide values
write_csv(norm_sel_coef, 'filtered_coefficient_max.csv')

# compute median score for each aa sequence
norm_sel_coef%>%
dplyr::group_by(aa_seq, pool, condition, family)%>%
dplyr::filter(timepoint_end == 'T2')%>%
dplyr::mutate(med_norm = median(norm_sel_coef, na.rm =TRUE), 
      med_sel = median(sel_coef, na.rm = TRUE))%>%
select(condition, aa_seq, family, timepoint_end, pool, assay, pos1, pos2, 
        pos, max1, max2, max, deg1, deg2, deg, double_mutant, med_norm, med_sel, peptide_max)%>%
ungroup()%>%
unique()-> med_score

med_score$peptide_max<-
  replace_na(med_score$peptide_max, FALSE)

# reformat the peptide library family names 
med_score[!grepl('PRM_0',  med_score$family), 'family']<- 'reference'
med_score$family<-
  gsub(pattern = 'PRM_0', '', med_score$family)

## Generate heatmap
# order by amino acid per chemical properties
# positive : R, H, K
# negative : D, E
# Polar : S, T, N, Q 
# hydrophobic small: A, V, I, L, 
# hydrophobic large : M, F, Y, W
# special : C, G, P

aa <- 
  c('R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'A', 'V', 'I','L', 'M', 'F', 'Y', 'W', 'C', 'G', 'P', '*' )

med_score$deg1<- 
    factor(med_score$deg1, levels = aa)
med_score$deg2<- 
    factor(med_score$deg2, levels = aa)
med_score$deg<- 
    factor(med_score$deg, levels = aa)

# loop through the pools to generate the propers heatmap
for(i in unique(med_score$pool)) {
  
  # select the pool data
  df <- med_score%>%
    subset(pool == i)
  
  # set pool specific variables
  pool <- i
  assay <- unique(df$assay)

  # set color themes for the heatmap according to the assay
  if(assay == 'Avail'){
    c_viridis = 'E'
    c_point = 'red'
    low_point='red'
   
  }else if (assay == 'PPI'){
    c_viridis = 'F'
    c_point = 'black'
    low_point = 'white'
    
  }

  # select data for the single mutant
  df%>%
    subset(!double_mutant | peptide_max)-> df_s
  

  ## Heatmap for single mutants
  if(nrow(df_s)> 50){
    
    # set width parameter for saving the graph
    w = 3*2.5
    # set height parameter for saving the graph
    if(length(unique(df_s$condition))>1){
      h = 6
    }else {
      h=4
    }

  # if there is MTX condition for this pool 
  if(c('MTX') %in% unique(df_s$condition)){ 
  # generate heatmap for normalized data (only for MTX condition)  
      df_s%>%
      group_by(pos, family) %>%
      mutate(
      med_group = ifelse(med_norm >= 0.35, "high", "low")
        )%>%
      filter(n() > 2 & condition == 'MTX') %>%
      ungroup()->x
        
      y<-x[is.na(x$pos), ]
      
      y[y$deg1 != y$max1, 'deg']<-y[y$deg1 != y$max1, 'deg1']
      y[y$deg1 != y$max1, 'pos']<-y[y$deg1 != y$max1, 'pos1']
      y[y$deg1 != y$max1, 'max']<-y[y$deg1 != y$max1, 'max1']
      
      y[y$deg2 != y$max2, 'deg']<-y[y$deg2 != y$max2, 'deg2']
      y[y$deg2 != y$max2, 'pos']<-y[y$deg2 != y$max2, 'pos2']
      y[y$deg2 != y$max2, 'max']<-y[y$deg2 != y$max2, 'max2']
     
      m1<-y[is.na(y$deg),]
      m1$pos<-m1$pos1
      m1$deg<-m1$deg1
      m1$max<-m1$max1
      
      m2<-y[is.na(y$deg),]
      m2$pos<-m2$pos2
      m2$deg<-m2$deg2
      m2$max<-m2$max2
      
      x%>%
        bind_rows(y)%>%
        bind_rows(m1)%>%
        bind_rows(m2)%>%
        drop_na(deg)%>%
        drop_na(pos)%>%
        drop_na(max)->x

  
      ggplot(x)+
      facet_grid(cols = vars(family), #rows = vars(condition),
                scales ='free')+
      geom_tile(aes(x = as.factor(pos), y =deg, fill = med_norm))+
      geom_point(data = x[x$peptide_max, ], aes(x = as.factor(pos), y =max, color = med_group), shape = 4)+
      xlab('position')+
      ylab('amino acid')+
      scale_fill_viridis_c(
        option = c_viridis,
        na.value = 'white',
        name = paste0(assay, ' score'),
        limits = c(0, 1),
        oob = scales::oob_squish
      )+
      scale_color_manual(values = c("low" = low_point, "high" = c_point))+
      t+
      theme(legend.position = 'bottom')+
    guides(color = 'none')->p_s
    
    saveRDS(p_s, paste0('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/figures/filter_outlier/', pool, '_', assay, '_norm_sel_single_both.rds'))

    ggsave(paste0('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/figures/filter_outlier/', pool, '_', assay, '_norm_sel_single_both.png'),
           p_s, height = 4, width = w)
    
 }

  # make sure there is enough data for a heatmap
  p_s <- 
    df_s%>%
    mutate(
      med_group = ifelse(med_sel >= -0.5, "high", "low")
        )%>%
    group_by(pos, family) %>%
    filter(n() > 15) %>%
    ungroup()%>%
    drop_na(deg)

  if(nrow(p_s)>10){
    # generate heatmap for selection coefficient
    p_s<-
    ggplot(p_s)+
    facet_grid(cols = vars(family), rows = vars(condition),
               scales ='free')+
    geom_tile(aes(x = as.factor(pos), y =deg, fill = med_sel))+
    geom_point(data = p_s[p_s$peptide_max, ], aes(x = as.factor(pos), y =max, color = med_group), shape = 4)+
    xlab('position')+
    ylab('amino acid')+
    scale_fill_viridis_c(option = c_viridis, na.value = 'white',  limits =c(-1, 1), 
    name = 'selection coefficient', oob = scales::oob_squish)+
    scale_color_manual(values = c("low" = low_point, "high" = c_point)) +
    t+
    theme(legend.position = 'bottom')+
    guides(color = 'none')

  
  ggsave(paste0('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/figures/filter_outlier/', pool, '_', assay, '_sel_single_both.png'), p_s, height = h, width = w)
  }
  }

  ## Heatmap for double mutants

  # if the pool does not have double mutants
   if(sum(df$double_mutant, na.rm = TRUE)< 20){
    next
   } 
   # set width parameter for saving the graph
  if(pool %in% c('P1', 'P4', 'P5', 'P6')){
    w = 3*2.5
  } else if (pool == 'P2') {
     w = 7*2.5
  } else if (pool == 'P3') {
     w = 4*2.5
  }

# # generate heatmap for selection coefficient

  df%>%
  mutate(
      med_group = ifelse(med_sel >= -0.5, "high", "low")
        )%>%
  group_by(pos1, pos2, family, condition) %>%
   
   filter(n() > 15) %>%
   drop_na(deg1)->x
  
  p_d<- 
    ggplot(x)+
    facet_grid(cols = vars(family), rows = vars(condition),
               scales ='free')+
    geom_tile(aes(x = deg1, y = deg2, fill = med_sel))+
    geom_point(data = x[x$peptide_max, ], 
              aes(x = max1, y =max2, color =med_group), shape = 4)+
    xlab('deg. position 1')+
    ylab('deg. position 2')+
    scale_fill_viridis_c(option = c_viridis, na.value = 'white',  limits =c(-1, 1), oob = scales::oob_squish,
    name = 'selection coefficient')+
    scale_color_manual(values = c("low" = low_point, "high" = c_point))+
    t+
    theme(legend.position = 'bottom')+
    guides(color = 'none')

ggsave(paste0('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/figures/filter_outlier/', pool, '_', assay, '_sel_double_both.png'), p_d, height = h, width = w)

if(c('MTX') %in% unique(df_s$condition)){
# generate heatmap for normalized data (only for MTX condition)  

  df%>%
  drop_na(deg1)%>%
   mutate(
    med_group = ifelse(med_norm >= 0.35, "high", "low")
  )%>%
  group_by(pos1, family) %>%
    filter(n() > 2 & condition == 'MTX') ->x
  
 p_d<- 
    ggplot(x)+
    facet_grid(cols = vars(family),# rows = vars(condition),
               scales ='free')+
    geom_tile(aes(x = as.factor(deg1), y = as.factor(deg2), fill = med_norm))+
    geom_point(data = x[x$peptide_max, ], aes(x = as.factor(max1), y =as.factor(max2), color = med_group), shape = 4)+
    xlab('deg. position 1')+
    ylab('deg. position 2')+
    scale_fill_viridis_c(option = c_viridis, na.value = 'white',  limits =c(0, 1), oob = scales::oob_squish,
    name =  paste0(assay, ' score'))+
    scale_color_manual(values = c("low" = low_point, "high" = c_point)) +
    t+
    theme(legend.position = 'bottom')+
    guides(color = 'none')
saveRDS(p_d, paste0('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/figures/filter_outlier/', pool, '_', assay, '_norm_double_both.rds'))
ggsave(paste0('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/figures/filter_outlier/', pool, '_', assay, '_norm_double_both.png'), p_d, height = 4, width = w)
  }
}



# Distribution of effect

med_score%<>%
  mutate(pbd = case_when(
  pool == 'P4' ~ 363,
  pool == 'P5' ~ 366,
  pool == 'P6' ~ 385))


med_score%<>%
  dplyr::mutate(paired= family == pbd)

med_score%<>%
  mutate(name = case_when(
  pool == 'P4' ~ 'PBD 363',
  pool == 'P5' ~ 'PBD 366',
  pool == 'P6' ~ 'PBD 385'))

MTX_dist <- 
med_score%>%
  subset(deg !='reference' & condition == 'MTX' & pool %in% c('P4', 'P5', 'P6') & pbd %in% c('363', '366', '385'))%>%
  drop_na(pool, condition)%>%
ggplot()+
facet_grid(cols = vars(name), 
           rows = vars(condition), scales = 'free')+
  new_scale_color() +   
  geom_violin(aes(y = med_norm, x = family, fill = paired), 
              trim = FALSE, width = 0.8, draw_quantiles = c(0.25,0.75))+
  scale_fill_manual(values = c('grey85','transparent'))+
  new_scale_color() +   
  geom_point(data = med_score[med_score$peptide_max & med_score$condition =='MTX' & med_score$pool != 'P1', ], 
            aes(y =med_norm, x=family, color = family), size=1.5)+
  geom_hline(aes(yintercept=0), linetype = 'dashed')+
  geom_hline(aes(yintercept=1), linetype = 'dashed')+
  #stat_compare_means(aes(med_norm, family),)+
  stat_summary(aes(y = med_norm, x=family, color = library), fun=median, geom="point", size=1.5, color = 'black', shape = 4)+
  scale_color_manual(values =c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  ylab('PPI score')+
  xlab('peptide family')+
  ylim(-0.5, 1.2)+
  t+
  theme(legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())

saveRDS(MTX_dist, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/MTX_dist.rds')

DMSO_dist<-
med_score%>%
  subset(condition == 'DMSO' & family !='reference' & pool %in% c('P4', 'P5', 'P6'))%>%
 dplyr::group_by(family, pool)%>%
  filter(n()>20)%>%
ggplot()+
facet_grid(cols = vars(name), rows = vars(condition), scales = 'free')+
  geom_violin(aes(y = med_sel, x = family, fill=paired), 
              trim = FALSE, width = 0.8, draw_quantiles = c(0.25,0.75))+
  stat_summary(aes(y = med_sel, x=family, color = library), fun=median, geom="point", size=1.5, color = 'black', shape = 4)+
  scale_fill_manual(values = c('grey80','transparent'))+
  new_scale_color() +     
  geom_point(data = subset(med_score, peptide_max & condition == 'DMSO' &
                          (pool %in% c( 'P4', 'P5', 'P6'))), 
            aes(y=med_sel, x =family, color = family), size=1.5)+
  stat_summary(aes(y = med_sel, x =family), 
               fun=median, geom="point", size=1.5, color = 'black', shape = 4)+
  scale_color_manual(values =c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF")[5:7])+
  ylab('selection coefficient')+
  xlab('peptide family')+
  ylim(-0.2, 0.15)+
  t+
  theme(legend.position = 'none', strip.background.x= element_blank(), 
        strip.text.x = element_blank())

saveRDS(DMSO_dist, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/DMSO_dist.rds')

plot_grid(MTX_dist, DMSO_dist, nrow = 2, 
          rel_heights = c(1,1))->supp_fitness

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/supp_fitness.png', 
       height = 6, width = 8)


# save median normalized score for further analysis
write_csv(med_score, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/med_score_filtered.csv')


