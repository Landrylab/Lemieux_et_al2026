# AUTHOR : Pascale Lemieux
### Filter the replicate values ###
# import packages and functions
library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggsci)
library(gridExtra)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# import reference sequences
all_ref<-read_csv('all_reference.csv')

all_ref%>%
  subset(!duplicated(aa_seq))%>%
  select(family, aa_seq, pos1, pos2, pos, double_mutant, deg1, deg2, deg, max1, max2, max)->all_ref_aa

#import selection coefficient (output of 2_c_selection_analysis.R)
sel_coef<-read_csv('sel_coefficient_both_20.csv')

# add condition info for Pool 1 to 3
sel_coef[sel_coef$pool == 'P1', 'condition']<- 'MTX'
sel_coef[sel_coef$pool == 'P2', 'condition']<- 'DMSO'
sel_coef[sel_coef$pool == 'P3', 'condition']<- 'DMSO'

sel_coef%>%
  dplyr::filter(condition=='DMSO')%>%
  ggplot()+
  facet_grid(rows = vars(timepoint_end), cols = vars(Pool), scales = 'free_y')+
  geom_density(aes(x = sel_coef, color = replicate))+
  scale_color_manual(values = c('R1' = 'black', 'R2' = '#8c57a2ff', 'R3' = '#cf4e9cff'))+
  xlab('selection coefficient')+
  lims(x = c(-1.2, 1))+
  t+
  theme(legend.position = 'bottom')->s_coef
saveRDS(s_coef, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/sel_coef_distribution.rds')


# merge sel. coefficients with reference sequences information
sel_coef<-
  merge(sel_coef, 
  all_ref_aa, 
  all.y =FALSE, 
  all.x = TRUE) 


all_ref_aa[!grepl('PRM_',  all_ref_aa$family), 'deg']<-
  all_ref_aa[!grepl('PRM_',  all_ref_aa$family), 'family']

all_ref_aa[!grepl('PRM_',  all_ref_aa$family), 'family']<-'reference'

ref_seq <- 
  all_ref_aa[all_ref_aa$family == 'reference', ]

# extract reads info for the sequences used for normalization
sel_coef %>%
  filter( aa_seq %in% ref_seq$aa_seq)->ref_sel_coef_df


ref_sel_coef_df%>%
subset(timepoint_end == 'T2')%>%
ggplot()+
  facet_grid(rows= vars(Pool), cols= vars(condition), scales = 'free')+
  geom_point(aes(x = sel_coef, y = family, color = replicate))+
  scale_color_manual(values = c('R1' = 'black', 'R2' = '#8c57a2ff', 'R3' = '#cf4e9cff'))+
  lims(x = c(-1, 1))+
  xlab('selection coefficient')+
  ylab('reference peptide')+
  t+
  theme(legend.position = 'bottom')


ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/selection_coefficient_reference.png', width = 5, height = 10)

# compute the normalized sel. coefficient using 2 reference sequences 
# previously characterized with high and low signal in either the PPI or Avail assay 
norm_sel_coef <- 
  vector('list', length = 6)
names(norm_sel_coef) <-unique(sel_coef$pool)

for (i in unique(sel_coef$pool)) {
  df  <-  sel_coef%>%
                subset(pool==i)
  ref_df <- ref_sel_coef_df%>%
                subset(pool==i)
  
  norm_sel_coef[[i]] <- 
    get_norm_selcoef(df, ref_df)
  
}

## Check replicates of the sel. coefficient and the normalized signal before filtering
df_norm_sel_coef<-
  bind_rows(norm_sel_coef)

df_norm_sel_coef%>%
  dplyr::filter(pool %in% c('P1', 'P4', 'P5', 'P6') & condition =='MTX')%>%
  ggplot()+
  facet_grid(cols = vars(Pool), scales = 'free_y')+
  geom_density(aes(x = norm_sel_coef, color = replicate))+
  scale_color_manual(values = c('R1' = 'black', 'R2' = '#8c57a2ff', 'R3' = '#cf4e9cff'))+
  xlab('normalized s')+
  lims(x = c(-1.2, 1))+
  t+
  theme(legend.position = 'bottom')->norm_coef_distribution_bf

df_norm_sel_coef%>%
  group_by(pool, condition)%>%
  dplyr::summarise(n())

saveRDS(norm_coef_distribution_bf, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/norm_coef_distribution_bf.rds')

# loop through the 6 pools
for (i in unique(df_norm_sel_coef$pool)) {

# custom function for vizualisation of correlation between replicates
cor_plot(df_norm_sel_coef, p = i, 
        file_dir = '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/', 
        file_name = 'sel_repcheck')

}

# Filter outliers in MTX condition (only pool P1, P4, P5 and P6) 
# since their is no interesting signal in DMSO condition
wide_df<- 
  df_norm_sel_coef%>%
  filter(timepoint_end == 'T2' & pool %in% c('P4', 'P5', 'P6', 'P1') & condition == 'MTX')%>%
  select(replicate, condition, aa_seq, sequence, family, timepoint_end, pool, norm_sel_coef)%>%
  group_by(pool, aa_seq, condition)%>%
  pivot_wider(
              names_from = replicate, 
              values_from = norm_sel_coef, values_fn = unique)%>%
              ungroup()

filtered_data<-vector(mode = 'list', length = 4)
names(filtered_data)<-c('P4', 'P5', 'P6', 'P1')

# loop through the pools of interest for filtering
for(p_name in c('P1','P4', 'P5', 'P6')){

sub_pool<-
  wide_df%>%
  filter(pool == p_name)

# step 1 : get linear model from the comparison of the 3 replicates
lm_R1R2<-glm(R1 ~ R2, data = sub_pool)
lm_R1R3<-glm(R1 ~ R3, data = sub_pool)
lm_R2R3<-glm(R2 ~ R3, data = sub_pool)

# step 2 : get the parameters of the lm and their standard errors
coef_error<-
  bind_cols(
  comp = c('R1-R2', 
          'R1-R3', 
          'R2-R3'),
  coef = 
  c(summary(lm_R1R2)$coefficients[2, 1],
  summary(lm_R1R3)$coefficients[2, 1],
  summary(lm_R2R3)$coefficients[2, 1]), 
  error_coef = 
  c(summary(lm_R1R2)$coefficients[2, 2],
  summary(lm_R1R3)$coefficients[2, 2],
  summary(lm_R2R3)$coefficients[2, 2]),

  intercept = 
  c(summary(lm_R1R2)$coefficients[1, 1],
  summary(lm_R1R3)$coefficients[1, 1],
  summary(lm_R2R3)$coefficients[1, 1]), 
  error_inter = 
  c(summary(lm_R1R2)$coefficients[1, 2],
  summary(lm_R1R3)$coefficients[1, 2],
  summary(lm_R2R3)$coefficients[1, 2])
  )

# step 3 : set the tolerance treshold based on the intercept and slop difference of the lm
# set a max value of tolerance of 0.5
coef_error%<>%
  dplyr::mutate(max_coef = coef + error_coef, 
                min_coef = coef - error_coef,
                max_inter = intercept + error_inter,
                min_inter = intercept - error_inter
                )

tolerance = 
  (max(coef_error$max_coef) - min(coef_error$min_coef))+(max(coef_error$max_inter) - min(coef_error$min_inter))

# was at 0.4
if(tolerance > 0.4){
  tolerance = 0.4
}

# step 4 : compute the absolute difference between each replicate
sub_pool$abs_R1vR2<- abs(sub_pool$R1 - sub_pool$R2) 
sub_pool$abs_R1vR3<- abs(sub_pool$R1 - sub_pool$R3)
sub_pool$abs_R2vR3<- abs(sub_pool$R3 - sub_pool$R2)

# step 5 : filter the data based on the abs. diff. and the tolerance treshold
sub_pool$keep_rows_R1vR2<- abs(sub_pool$R1 - sub_pool$R2) <= tolerance  
sub_pool$keep_rows_R1vR3<- abs(sub_pool$R1 - sub_pool$R3) <= tolerance 
sub_pool$keep_rows_R2vR3<- abs(sub_pool$R3 - sub_pool$R2) <= tolerance

sub_pool$keep_rows_R1vR2<-
  replace_na(sub_pool$keep_rows_R1vR2, TRUE)

sub_pool$keep_rows_R1vR3<-
  replace_na(sub_pool$keep_rows_R1vR3, TRUE)

sub_pool$keep_rows_R2vR3<-
  replace_na(sub_pool$keep_rows_R2vR3, TRUE)

# Identify which data point are respecting the tolerance for the comparisons between
# the 3 replicates
sub_pool$sum_keep<-
  rowSums(sub_pool[, c('keep_rows_R1vR2', 'keep_rows_R1vR3', 'keep_rows_R2vR3')])

# Visualise one comparison between a pair of replicate and the identify which data point respect the tolerance
sub_pool%>%
ggplot()+
geom_point(aes(x = R1, y = R3, color = as.factor(sum_keep)))+
lims(x = c(-0.5,1), y = c(-0.5,1))+
t->tolerance_graph
# Save the comparison
saveRDS(tolerance_graph, paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/', p_name, '_tolerance.rds'))

# step 6 : remove the replicate values for which the absolute difference is above the tolerance treshold 
# for 2/3 comparisons
sub_pool%>%
filter(sum_keep<=2)%>% ## was 1
select(sequence, family, pool, abs_R1vR2, abs_R1vR3, abs_R2vR3)%>%
pivot_longer(cols = c('abs_R1vR2', 'abs_R1vR3', 'abs_R2vR3'), 
              names_to = 'comparison', names_prefix = 'abs_', 
              values_to = 'abs_diff')->lonely_rep

# If there is no abs. difference comparison that is above the tolerance treshold, 
# save the correlation plot (same as above)
if(nrow(lonely_rep) == 0){
  sub_pool%<>%
  select(condition, aa_seq, sequence, family, timepoint_end, pool, R1, R2, R3)%>%
  pivot_longer(cols = c('R1', 'R2', 'R3'), 
            names_to = 'replicate', 
            values_to = 'norm_sel_coef')
  
  cor_plot(df =sub_pool, p = p_name, just_norm = TRUE,
        file_dir='~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/', 
        file_name = 'repcheck_filter')

print(p_name)
filtered_data[[p_name]]<-sub_pool
next
}              

# Identify which replicates are involed in the abs. difference that are above the treshold
lonely_rep<-bind_cols(lonely_rep, 
matrix(unlist(strsplit(lonely_rep$comparison, 
        split = 'v')), ncol = 2, byrow = TRUE))

colnames(lonely_rep)[6:7]<-c('r1', 'r2')

lonely_rep%<>%
dplyr::filter(abs_diff>tolerance)%>%
dplyr::group_by(sequence)%>%
dplyr::summarise(
    replicate = names(which.max(table(c(r1, r2))))
  )

# The data points have to be removed from the replicate
lonely_rep$remove<-TRUE

sub_pool%<>%
select(condition, aa_seq, sequence, family, timepoint_end, pool, R1, R2, R3)%>%
pivot_longer(cols = c('R1', 'R2', 'R3'), 
            names_to = 'replicate', 
            values_to = 'norm_sel_coef')%>%
  left_join( lonely_rep, by = c('sequence', 'replicate'))

sub_pool$remove<-
  replace_na(sub_pool$remove, FALSE)

sub_pool%>%
  filter(!remove)->filter_pool

# save the correlation plot for the filtered data
cor_plot(df =filter_pool, p = p_name, just_norm = TRUE,
        file_dir='~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/filter_outlier/', 
        file_name = 'repcheck_filter')

print(p_name)
filtered_data[[p_name]]<-filter_pool

}

filtered_data<-
  bind_rows(filtered_data)

filtered_data%>%
  group_by(pool, condition)%>%
  dplyr::summarise(n())


# join the full dataset by keeping the filtered observation in MTX condition
df_norm_sel_coef%>%
as_tibble()%>%
select(!c('norm_sel_coef'))%>%
  filter(condition == 'MTX')%>%
  right_join(filtered_data)->MTX_norm 

# Add the DMSO (unfiltered) observations
df_norm_sel_coef%>%
filter(condition == 'DMSO')%>%
bind_rows(MTX_norm)->filtered_data


filtered_data%>%
  dplyr::filter(pool %in% c('P1', 'P4', 'P5', 'P6') & condition =='MTX')%>%
  ggplot()+
  facet_grid(cols = vars(Pool), scales = 'free_y')+
  geom_density(aes(x = norm_sel_coef, color = replicate))+
  scale_color_manual(values = c('R1' = 'black', 'R2' = '#8c57a2ff', 'R3' = '#cf4e9cff'))+
  xlab('normalized s')+
  lims(x = c(-1.2, 1))+
  t+
  theme(legend.position = 'bottom')->norm_s

saveRDS(norm_s, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/norm_coef_distribution.rds')


# save the normalized filtered sel coefficient
write_csv(filtered_data, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/filtered_coefficient.csv')



