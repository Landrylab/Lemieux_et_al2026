# AUTHOR : Pascale Lemieux
# mutation additivity analysis i.e. epistasis

# import packages
library(tidyverse)
library(ggpubr)
library(magrittr)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/')

# import data
ppi_data <- read_csv('ppi_signif_score.csv')


# compare sum of individual mutation with the double mutations

ppi_data%>%
    dplyr::filter(max1 == deg1) %>%
    dplyr::select(aa_seq, family, deg2, med_norm_sel_coef, binding)-> position1
  
ppi_data%>%
  dplyr::filter(max2 == deg2)%>%
  dplyr::select(aa_seq, family, deg1, med_norm_sel_coef, binding)-> position2

ppi_data%>%
  dplyr::filter(max_pep)%>%
  dplyr::select(aa_seq, family, max1, max2, med_norm_sel_coef, binding) -> reference

# remove stop codons
ppi_data%>%
  dplyr::filter(max1 != deg1 & max2 != deg2 & deg1 != '*' & deg2 != '*')%>%
  dplyr::select(aa_seq, family, deg1, deg2, med_norm_sel_coef, binding) -> double

left_join(double, 
          position1, 
          by = c('family', 'deg2'), 
          suffix = c('.double', '.single2'))->single_vs_double

left_join(single_vs_double, 
          position2, 
          by = c('family', 'deg1'), 
          suffix = c('.double', '.single1'), 
          relationship = "many-to-many")->single_vs_double

colnames(single_vs_double)[10:12] <- c("aa_seq.single1", "med_norm_sel_coef.single1", "binding.single1")


left_join(single_vs_double, 
          reference[, c(2,5)], 
          by = 'family')-> single_vs_double

colnames(single_vs_double)[13] <- 'reference'

single_vs_double%<>%
  dplyr::group_by(aa_seq.double, aa_seq.single1, aa_seq.single2, reference)%>%
  dplyr::mutate(s1 = med_norm_sel_coef.single1-reference, 
                s2 = med_norm_sel_coef.single2-reference, 
                double = med_norm_sel_coef.double-reference)

single_vs_double%<>%
  dplyr::mutate(exp_epistasis = sum(s1, s2))


single_vs_double%>%
  dplyr::select(aa_seq.double, family, s1, s2, double)%>%
  pivot_longer(values_to = 'diff_score', cols = c(s1, s2, double), names_to = 'mutant')%>%
  ungroup()%>%
  select(aa_seq.double, diff_score, mutant)%>%
  unique()->dist_diff

ggplot(dist_diff)+
  facet_wrap(vars(family))+
  geom_histogram(aes(diff_score, fill =mutant), alpha = 0.5)+
  geom_vline(xintercept = 0)+
  t

  
single_vs_double%<>%
  dplyr::mutate(dif_epistasis = double + abs(exp_epistasis))

single_vs_double$name <- 
  paste0('PBD ', single_vs_double$family)


Fig2F <- 
single_vs_double%>%
  dplyr::filter(name == 'PBD 366')%>%
ggplot()+
  facet_wrap(vars(name))+
  geom_abline(linetype = 'dashed')+
  geom_point(aes(double, exp_epistasis, color = med_norm_sel_coef.double), alpha = 0.6)+
  stat_smooth(aes(double, exp_epistasis), method = 'loess', color ='black', fullrange=TRUE)+
  scale_color_viridis_c(option = 'F', direction = 1, end = 0.95)+
  stat_cor(aes(double, exp_epistasis), cor.coef.name = 'rho', size = 4, label.sep = '\n')+
  xlab(expression(paste(Delta, '(PPI'['max'], ' - exp. PPI'['double mutant'], ')')))+
  ylab(expression(paste(Delta, '(PPI'['max'], ' - pred. PPI'['double mutant'], ')')))+
  lims(x = c(-0.7, 0.2), y = c(-0.7, 0.2))+
  t+
  theme(strip.text = element_text(color = 'black'), 
        legend.position = 'none')

saveRDS(Fig2F, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Fig2F.rds')

single_vs_double%>%
  dplyr::filter(name != 'PBD 366')%>%
  ggplot()+
  geom_abline(linetype = 'dashed')+
  facet_wrap(vars(name))+
  geom_point(aes(double, exp_epistasis, color = med_norm_sel_coef.double), alpha = 0.6)+
  scale_color_viridis_c(option = 'F', direction = 1, end = 0.95, breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1))+
  stat_smooth(aes(double, exp_epistasis), method = 'auto', color ='black', fullrange=F)+
  stat_cor(aes(double, exp_epistasis), cor.coef.name = 'rho', size = 4, label.sep = '\n')+
  xlab(expression(paste(Delta, '(PPI'['max'], ' - exp. PPI'['double mutant'], ')')))+
  ylab(expression(paste(Delta, '(PPI'['max'], ' - pred. PPI'['double mutant'], ')')))+
  
  lims(x = c(-1.5, 0.3), y = c(-1.5, 0.3))+
  t+
  theme(strip.text = element_text(color = 'black'), 
        legend.position = 'bottom')+
  guides(color = guide_colorbar(title = 'PPI score', theme = theme(
    legend.key.width  = unit(8, "lines"),
    legend.key.height = unit(1, "lines")
  )))-> FigSAdditivity



saveRDS(FigSAdditivity, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Supp_additivity.rds')
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA1/Supp_additivity.png', FigSAdditivity, 
       width = 8, height = 6)



