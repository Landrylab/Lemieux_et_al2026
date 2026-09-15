# AUTHOR : Pascale Lemieux
# compare score PPI and avail score

# import libraries and custom functions

library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggupset)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# import median score from July2025 experiment
med_score <- 
  read_csv('signif_score_welch20_filter.csv')

med_score_pca1<-
  read_csv('../PCA1/ppi_signif_score.csv')

# get avail score from the first screen Feb2025
med_score_feb <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/avail_signif_score.csv')  

# format both datasets
med_score_feb$pos1 <- 
  as.numeric(med_score_feb$pos1)
med_score_feb$pos2 <- 
  as.numeric(med_score_feb$pos2)

med_score$family <-
  as.character(med_score$family)

# combine availability data from both experiments
med_score%>%
dplyr::filter(assay == 'Avail' & condition == 'MTX')-> med_avail

colnames(med_score_feb)[c(11,12)]<- c('med_norm', 'side')

med_avail <- bind_rows(med_score_feb, med_avail)

# select PPI score from the July2025 experiment
# select the designed PBD and peptide family combination 
med_score%>%
  dplyr::filter(assay == 'PPI' & condition == 'MTX')%>%
  dplyr::mutate(paired = case_when(family == '363' & pool == 'P4' ~ TRUE,
                                   family == '366' & pool == 'P5' ~ TRUE, 
                                   family == '385' & pool == 'P6' ~ TRUE))%>%
  dplyr::mutate(paired = replace_na(paired, FALSE))-> med_ppi


colnames(med_score_pca1)[c(9, 11,12)]<-c('peptide_max', 'med_norm', 'side')

med_ppi<-
  bind_rows( med_score_pca1, med_ppi, .id = 'screen')

# associate availability score with PPI score for each peptide sequence
merge(med_ppi[, c('aa_seq','family', 'med_norm', 'pool', 'side', 'screen', 'paired')], 
      med_avail[, c('aa_seq','family', 'med_norm', 'side')], 
      by.x = c('aa_seq', 'family'),
      by.y = c('aa_seq', 'family'),
      suffixes = c('.ppi', '.avail'), 
      all = FALSE) -> avail_ppi




# plot count of availability only for stronger interaction of the design paired of peptide family and PBD

avail_ppi$stop<-grepl('*', avail_ppi$aa_seq, fixed = TRUE)


#sub_pool$stop<-grepl('*', sub_pool$aa_seq, fixed = TRUE)

avail_ppi%>%
  dplyr::mutate(side.avail = factor(side.avail, levels = c('weaker','no significant difference', 'stronger'), 
                                    labels = c('weaker','no significant \ndifference', 'stronger')))%>%
  dplyr::filter(side.ppi == 'stronger' & !stop)%>%
  dplyr::filter(family %in% c(363, 366, 385))%>%
  ggplot()+
  geom_bar(aes(x = family, fill = side.avail), position ='fill', color = 'black')+
  t+
  scale_fill_manual(values = c( '#373737ff','#5b79a5ff', '#72c6ffff'))+
  theme(legend.position = 'none')+
  ylab('frequency')->stronger_ppi


avail_ppi%>%
  dplyr::mutate(side.avail = factor(side.avail, levels = c('weaker','no significant difference', 'stronger'), 
                                            labels = c('weaker','no significant \ndifference', 'stronger')))%>%
  dplyr::filter(screen ==2)%>%
  dplyr::filter(side.ppi == 'stronger' & !stop) %>%
  dplyr::group_by(side.avail)%>%
  dplyr::summarise(n())



124/(124+79+12)
#57% of stronger ppi are stronger avail

2326/(2326+3530+240)
#38% are stronger avail in all peptides

avail_ppi%>%
  group_by(side.ppi)%>%
  dplyr::filter(screen==2)%>%
  dplyr::summarise(n())

228/(5052+228+1254)

avail_ppi%>%
  dplyr::mutate(side.avail = factor(side.avail, levels = c('weaker','no significant difference', 'stronger'), 
                                    labels = c('weaker','no significant \ndifference', 'stronger')))%>%
  dplyr::filter(family %in% c(363, 366, 385) & !stop & paired)%>% 
  ggplot()+
  geom_bar(aes(x = family, fill = side.avail), position ='fill', color = 'black')+
  t+
  scale_fill_manual(values = c( '#373737ff','#5b79a5ff', '#72c6ffff'))+
  theme(legend.position = 'none')+
  ylab('frequency')->all_ppi

saveRDS(stronger_ppi, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/stronger_ppi_vs_avail.rds')
saveRDS(all_ppi, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/all_ppi_vs_avail.rds')

library(cowplot)

plot_grid(all_ppi, stronger_ppi)
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/ppi_vs_avail.png')


avail_ppi%>%
  #dplyr::filter(side.ppi == 'stronger')%>%
  dplyr::group_by(family, side.avail)%>%
  dplyr::summarise(n())


avail_ppi%<>%
  dplyr::mutate(paired = case_when(family == '363' & pool == 'P4' ~ TRUE,
                                   family == '366' & pool == 'P5' ~ TRUE, 
                                   family == '385' & pool == 'P6' ~ TRUE, 
                                   family %in% c('152', '250', '299', '246')~ TRUE))%>%
  dplyr::mutate(paired = replace_na(paired, FALSE))



# get hydrophobicity and PWM scores for each peptide variants

PWM_score_2<-
read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/PWM_score.csv')%>%
  select(aa_seq, family, pool, rel_PWM_score)

PWM_score_1<-
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/PWM_scores.csv')%>%
  select(aa_seq, rel_PWM_score, PBD)


PWM_score_1$PBD <-as.numeric(gsub(pattern = 'PBD ', '', PWM_score_1$PBD))
PWM_score_1$family <-as.numeric(PWM_score_1$PBD)

library(magrittr)

PWM_score_2%<>%
  dplyr::mutate(PBD = case_when(pool =='P4' ~ 363, 
                                pool == 'P5' ~ 366, 
                                pool =='P6' ~ 385))

all_PWM<-
  bind_rows(PWM_score_2[, -3], PWM_score_1)

# compute hydrophobicity and stickiness for all peptide

# import stickiness score for each amino acid
sticky <- read_table('~/PL_projects/PL_papers/PPI_optimization_paper/data/sticky_score.txt')
sticky2022<- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/stickiness2022.csv')
sticky_IDP<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/stickiness_IDP.csv')[, -1]

hydrophobicity<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/hydrophobicity.csv')
colnames(hydrophobicity)[1]<-'AA'

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

Cornette<-
  comp_sticky%>%
  select('AA', 'Cornette')

colnames(Cornette)[2]<-'sticky_score'

# identify stop codons
stop <- 
  grepl(all_PWM$aa_seq, pattern = '*', fixed = T)

# modify the aa sequence based on the stop codon
aa_seq <- 
  strsplit(all_PWM$aa_seq, '*', fixed = T)
aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')

# compute the sticky score for each aa sequence 

for(j in c('sticky', 'sticky2022', 'sticky_IDP', 'Cornette')){
  
  for (i in 1:length(aa_seq)) {
    if (length(aa_seq[[i]]) == 0){
      all_PWM[i, j] <- 0
      next
    }
    all_PWM[i, j] <- get_stickyness(aa_seq[[i]], get(j))
  }
}  


# see with lm which parameter explains more the PPI scores

mutatex<-
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/model_docking/ddG_interface.csv')%>%
  dplyr::filter(chain == 'B')%>%
  dplyr::mutate(PBD = as.character(PBD))%>%
  dplyr::group_by(PBD, chain, position, AA_mut, aa_ref)%>%
  dplyr::summarise(med_ddG = median(avg, na.rm =T))


med_score%>%
  dplyr::select(pool, aa_seq, family, deg, pos)%>%
  right_join(mutatex, join_by( deg==AA_mut,  family == PBD, pos ==position))%>%
  dplyr::select(aa_seq, med_ddG, family)%>%
  dplyr::mutate(family = as.numeric(family))%>%
  unique()->ddG_sequence



avail_ppi%>%
  dplyr::mutate(family = as.numeric(family))%>%
  dplyr::mutate(PBD = case_when(pool =='P4' ~ 363, 
                                pool == 'P5' ~ 366, 
                                pool =='P6' ~ 385, 
                                is.na(pool) ~ family))%>%
  
  left_join(all_PWM)%>%
  dplyr::filter(paired)%>%
  left_join(ddG_sequence)%>%
  unique()->test_model

med_ppi%>%
  dplyr::filter(peptide_max)%>%
  select(aa_seq)%>%
  unique()->ref_Seq


test_model$reference_seq<-
  test_model$aa_seq %in% unlist(ref_Seq)

write_csv2(test_model, '~/PL_projects/bootcamp_data.csv')


test_model%>%
  dplyr::filter(PBD == 385)->sub_data

PWM_model<-
  lm(med_norm.ppi ~ rel_PWM_score , data = sub_data)

summary(PWM_model)



ggplot(sub_data, aes(
  x = rel_PWM_score,
  y = med_norm.avail
)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_classic()




# select peptide to test with WesternBlot in the validation peptides already available
peptide_val<-read_csv('~/PL_projects/sequencing/Screen_Specificity_Fitness_July2025/peptide_validation_max.csv')

peptide_val$family<-
as.character(peptide_val$family)

left_join(peptide_val, 
          avail_ppi, 
          join_by(aa_seq, pool, family))->validation_avail_ppi


validation_avail_ppi%>%
  dplyr::group_by(family, pool)%>%
  filter(family == 385 )%>%
  #slice_max(order_by =med_norm.ppi, n=7)%>%
  ungroup()%>%
  select(aa_seq, family, med_norm.ppi, med_norm.avail, side.avail, side.ppi)


validation_avail_ppi%>%
  dplyr::group_by(family, pool)%>%
  slice_max(order_by =med_norm.avail, n=3)%>%
  ungroup()%>%
  select(aa_seq, family, med_norm.ppi, med_norm.avail)


colnames(sub_pool)[c(5,7)]<- c('ppi', 'avail')

# create a column to combine the results of stat test of avail and PPI
comp <- vector(mode='list', length = nrow(sub_pool))
for(i in 1:nrow(sub_pool)){
  
    x <- unlist(sub_pool[i, c('ppi', 'avail')])
    comp[[i]] <- str_c( x, names(x), sep = ' ')
  
}

sub_pool$comp <- comp


# Fisher's exact test
# keep only sequences with the data for the 2 assays and remove stop codons
comp <- 
  sub_pool[!grepl('*', sub_pool$aa_seq, fixed = TRUE), ]

# Initialize variables
m <- sum(table(comp$ppi)[c(2)])      # stronger PPI
n <- sum(table(comp$ppi)[c(1,3)])        # weaker/similar PPI
k <-  sum(table(comp$avail)[c(2)])       # stronger avail
x <-  1:sum(table(comp[, c('ppi', 'avail')])[c(2), c(2)]) # stronger PPI and stronger avail

# Use the dhyper built-in function for hypergeometric density
probabilities <- dhyper(x, m, n, k, log = FALSE)
probabilities

fisher.test(rbind(c(m,n),c(k,39)))

data <- data.frame( x = x, y = probabilities )
ggplot(data)+
  geom_density(stat="identity", aes(x, y))+
  geom_vline(xintercept = 39)

# Association is different from the null hypothesis (which is that both variable are independent)

### Do the same but with fitness assays vs PPI

# select PPI score from the July2025 experiment
med_score%>%
  filter(assay == 'PPI' & condition == 'MTX') -> med_ppi

# select fitness score from the July2025 experiment
med_score%>%
  filter(assay == 'PPI' & condition == 'DMSO') -> med_fitness


# associate availability score with PPI score for each peptide sequence
merge(med_ppi[, c('aa_seq','family', 'med_norm', 'pool', 'side')], 
      med_fitness[, c('aa_seq','family', 'med_sel', 'pool','side')], 
      by.x = c('aa_seq', 'family', 'pool'),
      by.y = c('aa_seq', 'family', 'pool'),
      suffixes = c('.ppi', '.fitness'), 
      all = FALSE) -> ppi_fitness


ppi_fitness%>%
  #filter(side.ppi == 'stronger')%>% 
  filter((pool == 'P4' & family == 363)|
           (pool == 'P5' & family == 366)|
           (pool == 'P6' & family == 385))%>%
  ggplot()+
  facet_grid(cols = vars(side.ppi), scales = 'free')+
  geom_bar(aes(x = side.fitness, fill = side.fitness), color = 'black')+
  scale_x_discrete(limits = c('stronger', 'no significant difference', 'weaker'), 
                   labels = c('stronger', 'no significant \ndifference', 'weaker'))+
  t+
  scale_fill_manual(values=c( 'white','#92D74D','#3D3576'))+
  theme(legend.position = 'none')

colnames(ppi_fitness)[c(5,7)]<- c('ppi', 'fitness')

# create a column to combine the results of stat test of avail and PPI
comp <- vector(mode='list', length = nrow(ppi_fitness))

for(i in 1:nrow(ppi_fitness)){
  
  x <- unlist(ppi_fitness[i, c('ppi', 'fitness')])
  comp[[i]] <- str_c( x, names(x), sep = ' ')
  
}

ppi_fitness$comp <- comp

# Fisher's exact test
# keep only sequences with the data for the 2 assays and remove stop codons
comp <- 
  ppi_fitness[!grepl('*', ppi_fitness$aa_seq, fixed = TRUE), ]

# Initialize variables
m <- sum(table(comp$ppi)['stronger'])      # stronger PPI
n <- sum(table(comp$ppi)[c('weaker', 'no significant difference')])        # weaker/similar PPI
k <-  sum(table(comp$fitness)['stronger'])       # stronger fitness
x <-  1:sum(table(comp[, c('ppi', 'fitness')])['stronger', 'stronger']) # stronger PPI and stronger fitness

# Use the dhyper built-in function for hypergeometric density
probabilities <- dhyper(x, m, n, k, log = FALSE)
probabilities

fisher.test(rbind(c(m,n),c(k,44)))

data <- data.frame( x = x, y = probabilities )
ggplot(data)+
  geom_density(stat="identity", aes(x, y))+
  geom_vline(xintercept = 44)

# there is an association between fitness and ppi strength also

