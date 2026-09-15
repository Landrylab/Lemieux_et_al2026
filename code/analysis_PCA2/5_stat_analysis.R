# AUTHOR : Pascale Lemieux
# Statistical analysis of differences in PCA signal 

# import libraries and custom function
library(tidyverse)
library(rstatix)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# import selection coefficient
sel_coef <- 
  read_csv('filtered_coefficient_max.csv')

# identify reference sequences
sel_coef[!grepl('PRM_0', sel_coef$family), 'family']<-'reference'
sel_coef$family<-
  gsub('PRM_0', '', sel_coef$family)

# based on phage display
sel_coef$peptide_max<-
  replace_na(sel_coef$peptide_max, FALSE) 

# format peptide family
sel_coef$family<-
  as.factor(sel_coef$family)

# check for max peptide in each family and verify statistical power 
sel_coef%>%
subset(family !='reference' )%>%#& condition == 'MTX')%>%#& pool %in% c('P1', 'P4', 'P5', 'P6'))%>%
select(pool, aa_seq, norm_sel_coef,peptide_max, family, sel_coef, condition, replicate)%>%
unique()-> unique_seq

unique_seq%>%
  filter((pool != 'P5' | replicate != 'R3'))->test_withoutP5R3

unique_seq%<>%
  dplyr::mutate(experiment = case_when(pool == 'P1' ~ 'Availability', 
                                       pool =='P4' ~ 'PPI ~ PBD 363', 
                                       pool =='P5' ~ 'PPI ~ PBD 366', 
                                       pool =='P6' ~ 'PPI ~ PBD 385'))


unique_seq%>%
#test_withoutP5R3%>%
filter(family %in% c(363, 366, 385) & condition =='MTX')%>%
ggplot()+
facet_wrap(~experiment, scales = 'free_x')+
  geom_jitter(aes(x = family, y = norm_sel_coef, color = peptide_max, alpha = peptide_max), shape = 1)+
  scale_color_manual(values =c('grey', 'red'), 
                     labels =  c('peptide variant', expression(paste('PWM'[max]))))+
  scale_alpha_manual(values = c(0.15,1))+
  ylab('normalized s')+
  ylim(-1, 1.2)+
  xlab('peptide family')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())+
  guides(alpha = 'none')->PWM_ref

saveRDS(PWM_ref, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/PWM_max_vs_variant.rds')

unique_seq%>%
#test_withoutP5R3%>%
  filter(condition =='DMSO')%>%
  ggplot()+
  facet_wrap(~pool, scales = 'free_x')+
  geom_jitter(aes(x = family, y = sel_coef, color = peptide_max, alpha = peptide_max), shape = 1)+
  scale_color_manual(values =c('grey', 'red'))+
  scale_alpha_manual(values = c(0.1,1))+
  ylab('sel coef')+
  #ylim(-1, 1.2)+
  xlab('peptide family')+
  t+
  theme(legend.position = 'bottom')

#ggsave('./figures/max_pep_distribution.png')

# identify the peptide family to compare
vec_family <- unique(sel_coef$family)[unique(sel_coef$family) != 'reference']


sel_coef%>%
  filter((pool != 'P5' | replicate != 'R3'))->test_withoutP5R3


#sel_coef<-test_withoutP5R3
# create a list with the relevant data for the statistical tests
list_comp <- list(
  list('pool' ='P4', 'df' = sel_coef[sel_coef$pool == 'P4' & sel_coef$condition == 'MTX', ], name_cond = 'P4_MTX', metric = 'norm_sel_coef', rel_family = c('363', '366', '385')),
  list('pool' ='P5', 'df' = sel_coef[sel_coef$pool == 'P5' & sel_coef$condition == 'MTX', ],  name_cond = 'P5_MTX', metric = 'norm_sel_coef', rel_family = c('363', '366', '385')),
  list('pool' ='P6', 'df' = sel_coef[sel_coef$pool == 'P6' & sel_coef$condition == 'MTX', ],  name_cond = 'P6_MTX', metric = 'norm_sel_coef', rel_family = c('363', '366', '385')), 
  list('pool' = 'P1', 'df' = sel_coef[sel_coef$pool == 'P1' & sel_coef$condition == 'MTX', ],  name_cond = 'P1_MTX', metric = 'norm_sel_coef', rel_family = c('363', '366', '385')), 
  list('pool' = 'P4', 'df' = sel_coef[sel_coef$pool == 'P4' & sel_coef$condition == 'DMSO', ],  name_cond = 'P4_DMSO', metric = 'sel_coef', rel_family = c('363', '366', '385')),
  list('pool' = 'P5', 'df' = sel_coef[sel_coef$pool == 'P5' & sel_coef$condition == 'DMSO', ],  name_cond = 'P5_DMSO', metric = 'sel_coef', rel_family = c('363', '366', '385')),
  list('pool' = 'P6', 'df' = sel_coef[sel_coef$pool == 'P6' & sel_coef$condition == 'DMSO', ],  name_cond = 'P6_DMSO', metric = 'sel_coef', rel_family = c('363', '366', '385')), 
  list('pool' = 'P2', 'df' = sel_coef[sel_coef$pool == 'P2' & sel_coef$condition == 'DMSO', ],  name_cond = 'P2_DMSO', metric = 'sel_coef', rel_family = c('152', '246', '250', '299', '363', '366', '385')),
  list('pool' = 'P3', 'df' = sel_coef[sel_coef$pool == 'P3' & sel_coef$condition == 'DMSO', ],  name_cond = 'P3_DMSO', metric = 'sel_coef', rel_family = c('152', '246', '250', '299'))
)

# create object to store the results
signif_df <- list()

#loop through the datasets of each pool
for (i in 1:length(list_comp)) {

  # remove the reference sequences from the pool 
  test_df <- list_comp[[i]]$df
  test_df%>%
    filter(family != 'reference') -> test_df

  pool <- list_comp[[i]]$pool
  # select only the relevant peptide family 
  vec_family <- unlist(unique(test_df[test_df$family %in% list_comp[[i]]$rel_family, 'family']))
  
  # test for difference between aa sequences (synonymous dna sequences) between the 
  # reference peptide (max binder found by phage display) and all the other sequences
  # using welch t-test with a Benjamini-Hochberg correction (FDR)
  stat_comp_intra_family(test_df, vec_family, welch = TRUE, metric=list_comp[[i]]$metric) -> signif_df[[list_comp[[i]]$name_cond]]

}

# combine the significant comparison from each pool
bind_rows(signif_df, .id ='condition') -> signif_all
signif_all$signif <- TRUE
colnames(signif_all)[2] <- 'aa_seq'

cond<-data.frame(matrix(unlist(strsplit(signif_all$condition, split = '_')), ncol =2, byrow = TRUE))

colnames(cond)<-c('pool', 'condition')

signif_all<-bind_cols(signif_all[, -1], cond)


# join the results of the tests with the PCA dataset
signif_data <- 
left_join(sel_coef, 
      signif_all, 
      by = c('aa_seq', 'family', 'pool', 'condition'))

signif_data$signif <- 
  replace_na(signif_data$signif, FALSE)

signif_data$side <- 
  replace_na(signif_data$side, 'no significant difference')


# compute the median PCA signal for the selection coef and the normalized scores.
signif_data%<>%
group_by(aa_seq, pool, condition)%>%
  dplyr::mutate(med_norm = median(norm_sel_coef, na.rm = TRUE),
         med_sel = median(sel_coef, na.rm=TRUE))%>%
  dplyr::select(pool, aa_seq,family, peptide_max, signif, med_norm, med_sel, deg1, deg2, pos1, pos2, max1, max2, deg, pos, max, assay, condition, side)%>%
  dplyr::distinct()

# save the results of the statistical tests 
write_csv(signif_data, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/signif_score_welch20_filter.csv')

## once test done once start here
# import data
signif_data <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/signif_score_welch20_filter.csv')

# format the classification of the stat test result and the peptide family
signif_data$side<- 
gsub(pattern = 'no significant difference', replacement = 'no significant\ndifference', signif_data$side)
signif_data$family<-
  as.factor(signif_data$family)

signif_data%<>%
  dplyr::mutate(experiment = case_when(pool == 'P1' ~ 'Availability', 
                                       pool =='P4' ~ 'PPI ~ PBD 363', 
                                       pool =='P5' ~ 'PPI ~ PBD 366', 
                                       pool =='P6' ~ 'PPI ~ PBD 385'))



# plot the significant differences vs the reference peptide
signif_data%>%
filter(family %in% c('363', '366', '385') & condition =='MTX')%>%
ggplot()+
  facet_wrap(~pool, scales = 'free_x')+
  geom_jitter(aes(x = family, y = med_norm, color = side, alpha = side), shape=1)+
  scale_color_manual(values = c( 'grey', '#4e67c8ff' , '#373737'))+
  scale_alpha_manual(values = c(0.3,1,0.2))+
  ylim(-1, 1.2)+
  ylab('PCA score')+
  xlab('PBD family')+
  theme_classic2()+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())+
  guides(color = guide_legend(nrow = 1))->signif_specificity

saveRDS(signif_specificity, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/stat_sel_coef_20_filter.rds')

# count the occurences of each test result category
signif_data%>%
  dplyr::group_by(family, pool, condition)%>%
  dplyr::count(side)%>%
  drop_na() -> count_binding


count_binding$side<-
factor(count_binding$side, levels = c('weaker', 'no significant\ndifference', 'stronger'), 
       labels = c('weaker', 'no significant\ndifference', 'stronger'))

## plot the count of significant differences vs the reference peptide
count_binding%>%
  dplyr::filter(family %in% c('363', '366', '385') & pool %in% c('P1', 'P4', 'P5', 'P6') & condition == 'MTX' )%>%
ggplot()+
  facet_wrap(~pool, scales = 'free_x')+
  geom_col(aes(x = family, y = n, fill = side), 
           color = 'black',  linewidth = 0.6, width = 0.5,  position = position_dodge2(width = 0.8, preserve = "single"))+
  scale_fill_manual(values = c( '#373737', 'grey', '#4e67c8ff'))+
  scale_y_continuous(transform = 'log2')+
  labs(x = 'peptide family', y = 'count')+
  t+
  theme(legend.position = 'bottom', legend.title = element_blank())->count_signif_plot

saveRDS(count_signif_plot, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/count_signif.rds')

# Visualize the increased PPI peptide specificity

signif_data%>%
  dplyr::filter(condition == 'MTX' & assay == 'PPI')%>%
  dplyr::filter(side == 'stronger')%>%
  dplyr::filter((family == 363 & pool == 'P4')|
                (family == 366 & pool == 'P5')|
                (family == 385 & pool == 'P6'))%>%
  select(aa_seq, family, side)%>%
  unique()->stronger_ppi


signif_data%>%
  dplyr::filter(condition == 'MTX' & assay == 'PPI')%>%
  dplyr::filter(side == 'stronger')%>%
  dplyr::filter((family == 363 & pool != 'P4')|
                  (family == 366 & pool != 'P5')|
                  (family == 385 & pool != 'P6'))%>%
  select(aa_seq, family, side)%>%
  unique()->stronger_unspe_ppi


signif_data%>%
  dplyr::filter(condition == 'MTX' & assay == 'PPI')%>%
  dplyr::filter(aa_seq %in% stronger_ppi$aa_seq)%>%
  select(pool, aa_seq, side, family)%>%
  unique()->comp_stronger

comp_stronger%<>%
  dplyr::mutate(PBD = case_when(
    pool=='P6'~385, 
    pool=='P5'~366, 
    pool=='P4'~363, 
  ))

comp_stronger$PBD<-
factor(comp_stronger$PBD, 
       levels = c(363, 366, 385), 
       labels = c('363', '366', '385'))


comp_stronger %<>%
  select(aa_seq, family, PBD, side)%>%
  pivot_wider(names_from = side, values_from = PBD, 
              values_fn = list)

comp_stronger%>%
  #select(stronger, `no significant\ndifference`, weaker)%>%
  dplyr::group_by(family, stronger)%>%
  dplyr::summarise(n())->sub

sub$spe<-
  unlist(lapply(sub$stronger, length))

table(sub[, c(3:4)])

8
(9+13+32+8)

library(ggupset)


FigStronger_pep<-
ggplot(comp_stronger, 
       aes(x = stronger)) +
  geom_bar(aes(fill = family), color ='black')+
  scale_fill_manual(values = c( "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  scale_x_upset(reverse = T, order_by = 'degree', 
                sets = c('363', '366', '385'))+
  xlab('stronger PPI')+
  t+
  theme(legend.position = 'bottom')+
  guides(fill = guide_legend(title = 'peptide family'))+
  theme_combmatrix(combmatrix.label.text = element_text(color = 'black'))
 
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/comp_stronger_pep.png')
saveRDS(FigStronger_pep, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/comp_stronger_pep.rds')

table(unlist(comp_stronger[,3]))

# check aa properties that lead to specific ppi increase

comp_stronger%>%
  left_join(signif_data, join_by(aa_seq, family))->check_properties


hydrophobicity<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/hydrophobicity.csv')
colnames(hydrophobicity)[1]<-'AA'

hydrophobicity_c<-hydrophobicity%>%
  select(AA, Cornette)

colnames(hydrophobicity_c)[2]<-'sticky_score'

aa_seq <- 
  strsplit(comp_stronger$aa_seq, '*', fixed = T)
aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')


for (i in 1:length(aa_seq)) {
  if (length(aa_seq[[i]]) == 0){
    comp_stronger[i, 'Cornette'] <- 0
    next
  }
  comp_stronger[i, 'Cornette'] <- get_stickyness(aa_seq[[i]], hydrophobicity_c)
}




aa_prop<-tibble::tibble(aa  = c( 'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q',  'A', 'V', 'I', 'L', 'M', 'F', 'Y','W',  'C', 'G', 'P'))

aa_prop%<>%
  dplyr::mutate(prop = case_when(aa %in% c('R', 'H', 'K', 'D', 'E') ~ 'charged',
                                 aa %in% c('S', 'T', 'N', 'Q') ~ 'polar', 
                                 aa %in% c( 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W') ~ 'hydrophobic', 
                                 aa %in% c( 'C', 'G', 'P') ~ 'special'), 
                size = case_when(aa %in% c('G', 'A', 'S') ~ 'small', 
                                 aa %in% c('V', 'P', 'C', 'T', 'D', 'N') ~ 'medium', 
                                 aa %in% c('F', 'Y', 'H', 'L', 'M', 'Q', 'W', 'I', 'R', 'E', 'K') ~ 'large'))


check_properties%<>%
  filter(assay == 'PPI', condition == 'MTX')%>%
  dplyr::mutate(double = case_when(deg1 == max1 | deg2 == max2 ~ FALSE, 
                                   deg1 !=max1 & deg2 != max2 ~ TRUE, 
                                   is.na(deg1) ~ FALSE))

check_properties%>%
  dplyr::select(aa_seq, family, stronger, double)%>%
  unique()->test_spe

test_spe$spe<-  
unlist(lapply(test_spe$stronger, FUN = length))==1



library(ggpattern)

test_spe%>%
  #dplyr::filter(spe)%>%
ggplot()+
  facet_grid(cols = vars(family))+
  geom_bar(aes(double, fill=spe), color = 'black')+
  scale_x_discrete(labels = c('single', 'double'))+
  scale_fill_manual(values =c('grey30','grey80'), labels = c('non-spe.', 'spe.'))+
  t+
  theme(axis.title.x = element_blank())+
  guides(fill = guide_legend(title = ''))-> single_double
 

check_properties%>%
  dplyr::filter(!is.na(deg))->single_lib

check_properties%>%
  dplyr::filter(is.na(deg))->double_lib

double_lib%<>%
  dplyr::mutate(double = case_when(deg1 == max1 | deg2 == max2 ~ FALSE, 
                                   deg1 !=max1 & deg2 != max2 ~ TRUE), 
                deg = case_when(deg1 == max1 & deg2 != max2 ~ deg2, 
                                deg1 !=max1 & deg2 == max2 ~ deg1), 
                max = case_when(deg1 == max1 & deg2 != max2 ~ max2, 
                                deg1 !=max1 & deg2 == max2 ~ max1), 
                pos = case_when(deg1 == max1 & deg2 != max2 ~ pos2, 
                                deg1 !=max1 & deg2 == max2 ~ pos1))

double_lib%>%
  dplyr::filter(!double)%>%
  dplyr::select(aa_seq, family, stronger, deg, max, pos)%>%
  unique()->single_in_double


single_lib%>%
  dplyr::select(aa_seq, family, stronger, deg, max, pos)%>%
  unique()%>%
  bind_rows(single_in_double)-> all_single

table(all_single$pos)

all_single%<>%
  left_join(aa_prop, join_by(deg == aa))%>%
  left_join(aa_prop, join_by(max == aa), suffix = c('.deg', '.max'))

all_single%<>%
  dplyr::mutate(same_size = size.deg == size.max, 
                same_prop = prop.deg == prop.max)


ggplot(all_single, aes(x = same_prop,
                       pattern = same_size))+
  scale_x_discrete(labels = c('diff.\n properties', 'same\n properties'))+
  geom_bar_pattern(
    color = 'black',
    fill ='white')+
  t+
  scale_pattern_manual(values = c('none', "crosshatch"))+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(), 
        axis.title.x = element_blank())#->sum_properties



aa_seq <- 
  strsplit(all_single$aa_seq, '*', fixed = T)
aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')


for (i in 1:length(aa_seq)) {
  if (length(aa_seq[[i]]) == 0){
    all_single[i, 'Cornette'] <- 0
    next
  }
  all_single[i, 'Cornette'] <- get_stickyness(aa_seq[[i]], hydrophobicity_c)
}


signif_data%>%
  dplyr::filter(condition == 'MTX' & assay == 'PPI')%>%
  dplyr::filter(peptide_max)%>%
  dplyr::filter((family == 363 & pool == 'P4')|
                  (family == 366 & pool == 'P5')|
                  (family == 385 & pool == 'P6'))%>%
  select(aa_seq, family, side)%>%
  unique()->ref_pep

aa_seq <- 
  strsplit(ref_pep$aa_seq, '*', fixed = T)
aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')

for (i in 1:length(aa_seq)) {
  if (length(aa_seq[[i]]) == 0){
    ref_pep[i, 'Cornette'] <- 0
    next
  }
  ref_pep[i, 'Cornette'] <- get_stickyness(aa_seq[[i]], hydrophobicity_c)
}


comp_stronger%<>%
  left_join(ref_pep[, c('family', 'Cornette')], join_by(family))

comp_stronger%<>%
  dplyr::group_by(family)%>%
  dplyr::mutate(rel.hydro = case_when(Cornette.x>=Cornette.y ~ 'more hydro', 
                                      Cornette.x<Cornette.y ~ 'less hydro'))%>%
  unique()


ggplot(comp_stronger)+
  geom_bar(aes(x = family, fill = rel.hydro), alpha = 0.5, position = position_dodge())



all_single%>%
  dplyr::select(aa_seq, family, stronger, deg, max, pos, prop.deg, prop.max)%>%
  pivot_longer(names_to = c('type_sub'), cols = c(prop.deg, prop.max), 
               values_to = c('type_prob'))->prob_sub


all_single%>%
  dplyr::select(aa_seq, family, stronger, deg, max, pos, size.deg, size.max)%>%
  pivot_longer(names_to = c('type_sub_2'), cols = c(size.deg, size.max), 
               values_to = c('type_size'))->size_sub

left_join(size_sub, prob_sub)%>%
  unique()


prob_sub$type_sub<-factor(prob_sub$type_sub, 
                          levels = c('prop.max', 'prop.deg'))





ggplot(comp_stronger, aes(x=rel.hydro, fill = family))+  
  geom_bar( color = 'black', 
           position = position_dodge2(preserve = "single"))+
  xlab('hydrophobicity')+
  t+
  scale_fill_manual(values = c( "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  scale_x_discrete(
                    labels = c('less', 'more'), 
                    limits = c('less hydro', 'more hydro'))+
  guides(pattern = guide_legend(title = '', position = 'bottom'))->aa_hydro


size_sub$type_sub_2<-factor(size_sub$type_sub_2, 
                          levels = c('size.max', 'size.deg'))

ggplot(size_sub, aes(x = type_size, fill = type_sub_2, pattern =type_sub_2))+  
  #geom_bar(, position = position_dodge())+
  geom_bar_pattern(
    color = 'black',
    fill ='white', 
    position = position_dodge2(preserve = 'single'))+
  xlab('aa size')+
  t+
  scale_pattern_manual(values = c( "crosshatch", 'none'), 
                       labels = c(expression(paste('PWM'[max])), 'stronger pep.'))+
  guides(pattern = guide_legend(title = '', position = 'bottom'))->size_prop



ggplot(all_single)+
  geom_bar(aes(x = pos, fill = family), color ='black')+
  scale_x_discrete(limits = c(1:10))+
  xlab('position')+
  scale_fill_manual(values = c( "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  t+
  theme(legend.title = element_blank())->pos_prop

l_pep<-
  get_legend(size_prop)

l_pos<-get_legend(pos_prop)

plot_grid(pos_prop+theme(legend.position = 'none'),
         size_prop+theme(legend.position = 'none'), 
          aa_prop+theme(legend.position = 'none'),
          nrow = 1, rel_widths = c(1,1, 1))-> bottom

plot_grid(single_double, sum_properties, ncols =2)

saveRDS(size_prop, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/size_prop.rds')
saveRDS(aa_prop, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/aa_prop.rds')
saveRDS(pos_prop, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/pos_prop.rds')
saveRDS(aa_hydro, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/aa_hydro.rds')

# compare with the first screen if the signif peptides are the same 
PPI_Feb2025<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ppi_signif_score.csv')

PPI_Feb2025%<>%
  filter(family %in% c('363', '366', '385'))%>%
  select(aa_seq, signif, binding, max_pep, family)


signif_data%>%
  filter((pool == 'P4' & family == 363)|
           (pool == 'P5' & family == 366)|
           (pool == 'P6' & family == 385))%>%
  filter(condition == 'MTX')%>%
  ungroup()%>%
select(family, aa_seq, condition, side)%>%
right_join(PPI_Feb2025, join_by(family, aa_seq))-> comp_signif_screen

colnames(comp_signif_screen)[c(4,6)]<-c('signif_spe', 'signif_ppi')

table(comp_signif_screen[, c('signif_spe', 'signif_ppi')])

comp_signif_screen$signif_spe<-gsub('\n', ' ', comp_signif_screen$signif_spe)

comp_signif_screen$signif_side<-
  comp_signif_screen$signif_spe == comp_signif_screen$signif_ppi

comp_signif_screen%>%
  dplyr::mutate(comp = str_c(signif_ppi, signif_spe, sep=' - '))->x

comp_signif_screen_g<-
ggplot(x)+
  stat_count(aes(x = fct_infreq(comp), fill = family))+
  scale_fill_manual(values = c("#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  xlab('ppi screen - specificity screen')+
  t+
  theme(axis.text.x = element_text(angle = 45, hjust =1, vjust =1))


saveRDS(comp_signif_screen_g, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/comp_signif_screens.rds')
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/comp_signif_screens.png')



sum(comp_signif_screen$signif_side)/nrow(comp_signif_screen)
# 58.5% of the comparisons are conserved with the treshold1 and with replicate 3 in Pool5
# 58% of the comparisons are conserved with the treshold1, without replicate 3 in Pool5
# 60.1% of the comparisons are conserved with the treshold2 and with replicate 3 in Pool5 -> keep the treshold at 0.4, but select that >=2 comparison respect the treshold
# 58% of the comparisons are conserved with the treshold 2, without replicate 3 in Pool5
# 58.8% of the comparisons are conserved with the treshold 2, filter 0.35, with replicate 3 in Pool5
# 58.6% of the comparisons are conserved with the treshold 1, filter 0.35, with replicate 3 in Pool5






