# AUTHOR : Pascale Lemieux
# import libraries and functions

library(tidyverse)
library(rstatix)
library(ggpubr)
library(Biobase)
library(ggnewscale)   

# import custom functions
source('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/code/functions.R')

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/data/')

# import significance scores
signif_data <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/Lemieux_et_al2026/data/PCA2/signif_score_welch20_filter.csv')

# select single mutants from single and double variants
signif_data%>%
    dplyr::filter((deg2 == max2) | (deg1 == max1 ) | !is.na(deg))%>%
    dplyr::filter(deg2 != '*' & deg1 != '*' | deg != '*')%>%
    dplyr::filter(condition == 'MTX')%>%
    dplyr::filter(pool %in% c('P4', 'P5', 'P6') & family %in% c(363, 366, 385))->single

unique(single$aa_seq)

# select the 2 sequences in the double mutant as background sequence with the highest significant score
signif_data%>%
dplyr::filter((deg2 != max2) | (deg1 != max1 ) & is.na(deg))%>%
dplyr::filter(deg2 != '*' & deg1 != '*')%>%
dplyr::filter(condition == 'MTX')%>%
dplyr::filter((family == '363' & pool == 'P4') |
           (family == '366' & pool == 'P5') |
           (family == '385' & pool == 'P6'))%>%
dplyr::group_by(pool, family)%>%
dplyr::filter(signif)%>%
slice_max(order_by = med_norm, n=2) ->top_double

# pivot longer the single mutation dataframe
single%>%
  dplyr::filter(is.na(pos))%>%
  dplyr::select(pool ,aa_seq, family, peptide_max, signif, med_norm, med_sel, assay, condition, side, deg1, deg2, max1, max2, pos1, pos2)%>%
 pivot_longer(
        cols = c(pos1, pos2),
        names_to = "pos_type",
        values_to = "pos_value"
    ) %>%
    pivot_longer(
        cols = c(max1, max2),
        names_to = "max_type",
        values_to = "max_value"
    ) ->x
    
# reformat pos and max information for selection
x$pos_type<-
    gsub(x$pos_type, pattern = 'pos', replacement = '') 
x$max_type<-
    gsub(x$max_type, pattern = 'max', replacement = '')

# filter the entries with the corresponding pos and max
x%>%
dplyr::filter(pos_type == max_type)%>%
    pivot_longer(
        cols = c(deg1, deg2),
        names_to = "deg_type",
        values_to = "deg_value"
    )->x

# reformat deg information for selection    
x$deg_type<-
    gsub(x$deg_type, pattern = 'deg', replacement = '')

# filter the variants encoding the max peptide (defined by phage display)
x%>%
  dplyr::filter(deg_type == max_type)%>%
  dplyr::filter(peptide_max)%>%
  dplyr::select(pool ,aa_seq, family, peptide_max, signif, med_norm, med_sel, assay, condition, side, deg_value, max_value, pos_value)%>%
unique()->max_pep

# select single mutants from the double variant library
#x%>%
#filter(deg_type == max_type)%>%
#filter(deg_value!=max_value)%>%
#select(pool ,aa_seq, family, peptide_max, signif, med_norm, med_sel, assay, condition, side, deg_value, max_value, pos_value)%>%
#unique()->single_from_double

# select the single mutants from single variant library
single%>%
  dplyr::filter(!is.na(pos))%>%
  dplyr::select(pool ,aa_seq, family, peptide_max, signif, med_norm, med_sel, assay, condition, side, deg, max, pos)->single_from_single

colnames(single_from_single)[11:13]<-c('deg_value', 'max_value', 'pos_value')

# combine max peptide information with single mutant info
all_single<-
  bind_rows(max_pep, single_from_single)

table(all_single$aa_seq)

# select 2 single mutants per pool and position with the highest med_norm
all_single%>%
    dplyr::group_by(pool, family, pos_value)%>%
    slice_max(order_by =med_norm, n =2)%>% 
  dplyr::mutate(
    top_label = if_else(row_number() == 1, "top1", 'top2')
  )-> top2_bypos

top2_bypos%<>%
  dplyr::filter((family == '363' & pool == 'P4') |
           (family == '366' & pool == 'P5') |
           (family == '385' & pool == 'P6'))

# compute the difference in score between the top 2 mutants per position and family
# select the 5 positions with the smallest difference between the top 2 mutants at this position
# then, select the top 2 positions with the strongest top1 score which is not the max peptide (according to phage display)
top2_bypos%>%
select(pool, aa_seq, family, top_label, med_norm, deg_value, pos_value, max_value)%>%
dplyr::group_by(family, pos_value)%>%
pivot_wider(values_from = c('med_norm', 'aa_seq', 'deg_value', 'max_value'), 
            names_from = top_label)%>%
dplyr::mutate(dif_top = med_norm_top1 - med_norm_top2)%>%
dplyr::ungroup()%>%
dplyr::group_by(family)%>%
  dplyr::filter(deg_value_top1 != max_value_top1)%>%
slice_min(order_by = dif_top, n=5)%>%
slice_max(order_by = med_norm_top1, n=2)->pos_choice

pos_choice%>%
  dplyr::select(family, pool, aa_seq_top1, deg_value_top1, max_value_top1, pos_value, med_norm_top1)%>%
  dplyr::mutate(top_label ='top1')->top1

pos_choice%>%
  dplyr::select(family, pool, aa_seq_top2, deg_value_top2, max_value_top2, pos_value,  med_norm_top2)%>%
  dplyr::mutate(top_label ='top2')->top2

colnames(top1)<-
  gsub('_top1', '', colnames(top1))

colnames(top2)<-
  gsub('_top2', '', colnames(top2))

bind_rows(top1, top2)->summary_top

write_csv(summary_top, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/summary_opt_single.csv')

pos_choice%>%
  dplyr::select(family, pool, aa_seq_top1, aa_seq_top2, pos_value)%>%
  pivot_longer(cols = c('aa_seq_top1', 'aa_seq_top2'))-> pos_choice

top2_bypos%>%
right_join(pos_choice)%>%
  dplyr::select(pool, family,pos_value, deg_value, top_label, med_norm)%>%
  unique()-> pos_choice

# select the aa_seq of the max variants from double library as the 2 background 
# to combine with the top 2 positions with single mutants
top_double%>%
select(aa_seq, family, med_norm)%>%
ungroup()%>%
unique()-> max_seq

# function to generate the sequence variants with all the combination of background sequences and single mutants

##  compute the PCA score while generating the PCA opt peptide
get_variants_max_peptide<-
function(pos_choice, max_seq) {

seq_variant <-vector(mode = 'character')
score_variant<-vector(mode ='numeric')

for(n in 1:nrow(max_seq)){
  
  f_name<-unlist(max_seq[n, 'family'])

df_pos<-pos_choice%>%
    filter(family == f_name)

unlist(max_seq[n, 'aa_seq'])->max_ref

df_pos%>%
  dplyr::select(deg_value, pos_value, med_norm)%>%
  ungroup()%>%
  pivot_wider(values_from = 'deg_value', 
            names_from = 'pos_value')->wide_pos

aa_pos<-
  expand.grid(wide_pos[, 4:5])%>%
  drop_na()

aa_pos<-
bind_cols(
'pos1' =  as.character(aa_pos[,1]),
'pos2' = as.character(aa_pos[,2]))

df_pos%>%
  ungroup()%>%
  dplyr::mutate(pos_id = case_when(pos_value == max(pos_value) ~ 'pos2', 
                                   pos_value == min(pos_value) ~ 'pos1'))%>%
  select(deg_value, med_norm, pos_id)->score_single


for(i in 1:nrow(aa_pos)){

split_pep<-
    unlist(strsplit(max_ref, split = ''))

pos<-unique(as.numeric(df_pos$pos_value))

split_pep[pos]<-aa_pos[i, ]

seq<-str_c(split_pep, collapse ='')

seq_variant<-c(seq_variant, seq)

}

for (k in 1:nrow(aa_pos)) {
  
  aa_pos[k, ]%>%
    pivot_longer(cols=c('pos1','pos2'))%>%
    ungroup()%>%
    left_join(score_single, join_by(name == pos_id, value == deg_value))->score_mutations
  
    score_variant_x<-unlist(max_seq[n, 'med_norm']+ sum(score_mutations$med_norm))
    
    score_variant<-c(score_variant, score_variant_x)
    
}

l<-length(seq_variant)    
print(l)
}

info_variant<-
  list(seq_variant, score_variant)

return(info_variant)
}

# get all variants with up to 4 positions mutated (from the max peptide defined by phage display)
x<-get_variants_max_peptide(pos_choice, max_seq)

x<-tibble('aa_seq'= x[[1]],  'PCA_score' = x[[2]])

f_variants<-bind_cols(fam = rep(unique(pos_choice$family), each=8), 
                aa_seq  = x)


# save aa sequences for validation
write_csv(f_variants, 'PCA_optimized_validation.csv')

write_csv(pos_choice, 'PCA2/max_single.csv')
write_csv(max_seq, 'PCA2/double_background.csv')


# visualize PCA optimization 
norm_sel_coef<-read_csv('PCA2/filtered_coefficient_max.csv')

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


med_score%>%
  dplyr::filter(double_mutant)%>%
  select(family, pos1, pos2)%>%
  unique()%>%
  pivot_longer(cols = c('pos1', 'pos2'), 
               names_to = 'id_pos', 
               values_to = 'pos')->double_pos

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

pbd<-
  tibble('pbd' = c(363, 366, 385), 'pool' = c('P4', 'P5', 'P6'))

# loop through the pools to generate the propers heatmap
for(i in c('P4', 'P5', 'P6')) {  
  # select the pool data
  df <- med_score%>%
    subset(pool == i)
  
  # set pool specific variables
  pool <- i
  assay <- unique(df$assay)
  
  # set color themes for the heatmap according to the assay
    c_viridis = 'F'
    c_point = 'black'
    low_point = 'white'
    
    
    
    

  # select data for the single mutant
  df%>%
    subset(!double_mutant | peptide_max)-> df_s
  
  
  ## Heatmap for single mutants

    # set width parameter for saving the graph
    w = 3*2.5
    # set height parameter for saving the graph
    if(length(unique(df_s$condition))>1){
      h = 6
    }else {
      h=4
    }
    
    # if there is MTX condition for this pool 

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
      
      

      pos_choice->single_opt
      
      top_double->double_background

       x->
        test
      
     test%<>%
       left_join(single_opt[,c(1:5)], 
                 join_by(deg == deg_value, pos == pos_value, family, pool))
      
    double_background%<>%
      select(deg1, deg2, pos1, pos2)%>%
      pivot_longer(cols = c(deg1, deg2), values_to = 'deg', names_to = 'deg_id')%>%
      pivot_longer(cols= c(pos1, pos2), values_to = 'pos', names_to = 'pos_id')%>%
      dplyr::filter((pos_id == 'pos1' & deg_id == 'deg1') |
                      (pos_id == 'pos2' & deg_id == 'deg2' ))
    
    
    
    double_background$background<-TRUE
    
    
     test%>%
      left_join(double_background[, c('deg', 'pos', 'background', 'family', 'pool')])->y
    
      
     ggplot(test)+
        facet_grid(cols = vars(family), #rows = vars(condition),
                   scales ='free')+
        geom_tile(aes(x = as.factor(pos), y =deg, fill = med_norm))+
        #geom_tile(data = subset(test, background), 
         #        aes(x = as.factor(pos), y =deg), color = 'black', fill = 'transparent', size = 0.8)+
        geom_tile(data = subset(test, !is.na(top_label)),
                  aes(x = as.factor(pos), y =deg, color = as.factor(top_label)), fill = 'transparent', size = 0.8)+
        scale_color_manual(values = c('#8fba74ff', '#7b00a6ff', 'transparent'))+
       guides(color = guide_legend(title = ''))+
        new_scale_color() +
        geom_point(data = x, aes(x = as.factor(pos), y =max, color = med_group), shape = 4)+
        scale_y_discrete(limits = aa)+
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
      

      
      saveRDS(p_s, paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/', pool, '_', assay, '_norm_sel_single_both.rds'))
      
      ggsave(paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/', pool, '_', assay, '_norm_sel_single_both.png'),
             p_s, height = 4, width = w)
    
      
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
      
      ggsave(paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/', pool, '_', assay, '_sel_double_both.png'), p_d, height = h, width = w)
      
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
          geom_tile(data = subset(x, aa_seq %in% top_double$aa_seq & family == unlist(pbd[pbd$pool==i, 'pbd'])),
                    aes(x = deg1, y =deg2), fill = 'transparent', size = 0.8, color = 'black')+
          scale_color_manual(values = c("low" = low_point, "high" = c_point)) +
          t+
          theme(legend.position = 'bottom')+
          guides(color = 'none')
        saveRDS(p_d, paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/', pool, '_', assay, '_norm_double_both.rds'))
        ggsave(paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/', pool, '_', assay, '_norm_double_both.png'), p_d, height = 4, width = w)
      }
      
        
}      

P4_double_norm<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4_PPI_norm_double_both.rds')
  
P5_double_norm<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P5_PPI_norm_double_both.rds')

P6_double_norm<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P6_PPI_norm_double_both.rds')

P4_single_norm<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P4_PPI_norm_sel_single_both.rds')

P5_single_norm<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P5_PPI_norm_sel_single_both.rds')

P6_single_norm<-
  readRDS('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/P6_PPI_norm_sel_single_both.rds')

P4<-
plot_grid(P4_double_norm+theme(legend.position = 'none'),
          P4_single_norm, nrow = 2, rel_heights = c(1, 1.2), 
          labels = c('A', 'B'), label_fontface = 'plain', label_size = 14)

P5<-
  plot_grid(P5_double_norm+theme(legend.position = 'none'),
            P5_single_norm, nrow = 2, rel_heights = c(1, 1.2), 
            labels = c('A', 'B'), label_fontface = 'plain', label_size = 14)

P6<-
  plot_grid(P6_double_norm+theme(legend.position = 'none'),
            P6_single_norm, nrow = 2, rel_heights = c(1, 1.2), 
            labels = c('A', 'B'), label_fontface = 'plain', label_size = 14)

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS10.png', P4, width=8, height = 7)

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS11.png', P5, width=8, height = 7)

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS12.png', P6, width=8, height = 7)

