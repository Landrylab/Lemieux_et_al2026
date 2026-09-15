# AUTHOR : Pascale Lemieux
# Analysis of stronger avail peptides

# import libraries

library(tidyverse)
library(rstatix)
library(ggpubr)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/')

# import selection coefficient
sel_coef <- 
  read_csv('norm_avail_sel_coef.csv')

sel_coef[!(grepl('PRM', sel_coef$family)), 'family']<-'reference'

sel_coef[(grepl('PRM', sel_coef$family)), 'family'] <- 
  gsub('PRM_0', '', unlist(sel_coef[(grepl('PRM', sel_coef$family)), 'family']))

# use T2 selection coefficient

T2_sel_coef <- 
  sel_coef[sel_coef$timepoint_end == 'T2', ]

T2_sel_coef$max_pep <- 
  c(T2_sel_coef$deg1 == T2_sel_coef$max1 & T2_sel_coef$deg2 == T2_sel_coef$max2)

T2_sel_coef$max_pep <- 
  replace_na(T2_sel_coef$max_pep, FALSE)


# check for max peptide in each family and verify statistic power between families
FigS4D <- 
ggplot(T2_sel_coef[T2_sel_coef$family !='reference', ])+
  geom_jitter(aes(x = family, y = norm_sel_coef), alpha = 0.15, shape = 1)+
  geom_jitter(data = T2_sel_coef[T2_sel_coef$max_pep & T2_sel_coef$family !='reference', ], 
              aes(x = family, y = norm_sel_coef, color = sequence))+
  scale_color_viridis_d(option = 'H')+
  ylab('normalized s')+
  xlab('peptide family')+
  t+theme(legend.position = 'none')

saveRDS(FigS4D, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4D.rds')

family_test <- 
  t_test(T2_sel_coef[T2_sel_coef$max_pep & T2_sel_coef$family !='reference', ], norm_sel_coef ~ family, p.adjust.method = 'BH')

my_breaks <- c(1, 1e-2 ,1e-6, 1e-12)

# look at the differences between peptide family
ggplot(family_test)+
  geom_tile(aes(group1, group2, fill = p.adj))+
  geom_text(data = family_test[family_test$p.adj.signif =='ns', ], 
            mapping = aes(group1, group2, label= p.adj.signif))+
  scale_fill_viridis_c(breaks = my_breaks, labels =my_breaks, trans = 'log2', option = 'F')+
  theme_classic2()+
  theme(axis.text = element_text(color = 'black'), 
        legend.position = c(0.7, 0.2), 
        legend.title = element_text(), 
        axis.title = element_blank())+
  guides(fill = guide_colorbar(limits = c(10, 1e-8), 
                               position = 'inside', 
                               title = 'p-value adjusted', 
                               direction = 'horizontal', 
                               title.position = 'top', 
                               theme = theme(legend.key.width = unit(12, 'lines'))))


vec_family <- unique(T2_sel_coef$family)[-5]
signif_all <- tibble()
for (family in vec_family){
  
  sub_sel_coef <- sel_coef[sel_coef$family == family, ]
 
   aa_to_compare <- 
    unique(sub_sel_coef$aa_seq)
  
  aa_ref <- 
    unlist(unique(sub_sel_coef[sub_sel_coef$max_pep, 'aa_seq' ]))
  
  aa_seq_stat <- 
    t_test(sub_sel_coef, norm_sel_coef ~ aa_seq, ref.group = aa_ref, p.adjust.method = 'BH')
  
 
  signif_diff<-  aa_seq_stat[aa_seq_stat$p.adj.signif != 'ns', 'group2']
  
  stronger <- signif_diff$group2 %in% unlist(aa_seq_stat[aa_seq_stat$statistic < 0 & aa_seq_stat$p.adj.signif != 'ns', 'group2'])
  weaker <- signif_diff$group2 %in% unlist(aa_seq_stat[aa_seq_stat$statistic > 0 & aa_seq_stat$p.adj.signif != 'ns', 'group2'])
  
  signif_diff$family <- family
  signif_diff[weaker, 'binding'] <- 'weaker'
  signif_diff[stronger, 'binding'] <- 'stronger'
  
  signif_all <- bind_rows(signif_all, signif_diff)
  
}
  
signif_all$signif <- TRUE

colnames(signif_all)[c(1,3)] <- c('aa_seq', 'availability')

signif_data <- 
left_join(T2_sel_coef, 
      signif_all, 
      by = c('aa_seq', 'family'))

signif_data$signif <- 
  replace_na(signif_data$signif, FALSE)

signif_data$availability <- 
  replace_na(signif_data$availability, 'no significant difference')

signif_data%>%
  dplyr::select(aa_seq, deg1, deg2, pos1, pos2, max1, max2, family, norm_sel_coef, max_pep, signif, availability)%>%
  dplyr::group_by(aa_seq)%>%
  dplyr::mutate(med_norm_sel_coef = median(norm_sel_coef))%>%
  dplyr::select(aa_seq, deg1, deg2, pos1, pos2, max1, max2, family, max_pep, signif, med_norm_sel_coef, availability)%>%
  dplyr::distinct() -> signif_data


FigS4H <- 
ggplot(signif_data[signif_data$family != 'reference', ])+
  geom_jitter(aes(x = family, y = med_norm_sel_coef, color = availability, alpha = availability), shape = 1)+
  scale_color_manual(values= c('grey70','#92D74D', '#3D3576'))+
  scale_alpha_manual(values = c(0.4,1,0.4))+
  ylab('Availability score')+
  xlab('PBD family')+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())

saveRDS(FigS4H, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4G.rds')

write_csv(signif_data, 'avail_signif_score.csv')




