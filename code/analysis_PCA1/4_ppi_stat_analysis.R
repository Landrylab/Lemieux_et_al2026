# AUTHOR : Pascale Lemieux
# Analysis of stronger binding peptides

# import libraries

library(tidyverse)
library(rstatix)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/')

# import selection coefficient
sel_coef <- 
  read_csv('norm_sel_coef_ppi.csv')


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
FigS4G <- 
ggplot(T2_sel_coef[T2_sel_coef$family !='reference', ])+
  geom_jitter(aes(x = family, y = norm_sel_coef), alpha = 0.15, shape = 1)+
  geom_jitter(data = T2_sel_coef[T2_sel_coef$max_pep & T2_sel_coef$family !='reference', ], 
              aes(x = family, y = norm_sel_coef, color = sequence))+
  scale_color_viridis_d(option = 'H')+
  ylab('normalized s')+
  xlab('peptide family')+
  t+
  theme(legend.position = 'none')

saveRDS(FigS4G, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4H.rds')

family_test <- 
  t_test(T2_sel_coef[T2_sel_coef$max_pep, ], norm_sel_coef ~ family, p.adjust.method = 'BH')

my_breaks <- c(1, 1e-2 ,1e-6, 1e-12)

# check differences between families
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

# test for difference between aa sequences (synonymous codons)
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

colnames(signif_all)[1] <- 'aa_seq'

signif_data <- 
left_join(T2_sel_coef, 
      signif_all, 
      by = c('aa_seq', 'family'))

signif_data$signif <- 
  replace_na(signif_data$signif, FALSE)

signif_data$binding <- 
  replace_na(signif_data$binding, 'no significant difference')

signif_data%>%
  dplyr::select(aa_seq, deg1, deg2, pos1, pos2, max1, max2, family, norm_sel_coef, max_pep, signif, binding)%>%
  dplyr::group_by(aa_seq)%>%
  dplyr::mutate(med_norm_sel_coef = median(norm_sel_coef))%>%
  dplyr::select(aa_seq, deg1, deg2, pos1, pos2, max1, max2, family, max_pep, signif, med_norm_sel_coef, binding)%>%
  dplyr::distinct() -> signif_data

write_csv(signif_data, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ppi_signif_score.csv')
signif_data <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/ppi_signif_score.csv')

signif_data$binding <- 
gsub(pattern = 'no significant difference', replacement = 'no significant\ndifference', signif_data$binding)

FigS4B <- 
ggplot(signif_data[signif_data$family != 'reference', ])+
  geom_jitter(aes(x = family, y = med_norm_sel_coef, color = binding, alpha = binding), shape=1)+
  scale_color_manual(values= c('grey', '#FE7F2D', '#05004E'))+
  scale_alpha_manual(values = c(0.4,1,0.4))+
  ylab('PPI score')+
  xlab('PBD family')+
  theme_classic2()+
  t+
  theme(legend.position = 'bottom', 
        legend.title = element_blank())+
  guides(color = guide_legend(nrow = 1))


saveRDS(FigS4B, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/FigS4I.rds')



signif_data$name <- 
  gsub(signif_data$family, pattern = 'PRM_0', replacement= '')

signif_data%>%
  dplyr::group_by(name)%>%
  dplyr::count(binding) -> count_binding

count_binding$binding <- 
factor(count_binding$binding, levels = c('weaker', 'no significant\ndifference', 'stronger'), 
       labels = c('weaker', 'no significant\n difference', 'stronger'))

Fig2D <- 
ggplot(count_binding[count_binding$name != 'reference', ])+
  geom_col(aes(x = name, y = n, fill = binding), 
           color = 'black',  linewidth = 0.6, width = 0.5,  position = position_dodge2(width = 0.8, preserve = "single"))+
  scale_fill_manual(values = c( '#05004E', '#F5FBE6','#FE7F2D'))+
  scale_y_continuous(transform = 'log2')+
  labs(x = 'peptide family', y = 'count')+
  t+
  theme(legend.position = 'bottom')

saveRDS(Fig2D, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/fig2D.rds')

