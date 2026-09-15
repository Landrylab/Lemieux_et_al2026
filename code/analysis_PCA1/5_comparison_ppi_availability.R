# AUTHOR : Pascale Lemieux
# compare PPI score and availability score

# import libraries

library(tidyverse)
library(rstatix)
library(ggupset)
library(ggpubr)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/')

med_avail <- 
  read_csv('avail_signif_score.csv')

med_avail$pos1 <- 
  as.numeric(med_avail$pos1)
med_avail$pos2 <- 
  as.numeric(med_avail$pos2)

med_ppi <- 
  read_csv('ppi_signif_score.csv')


comp_signif <- 
full_join(med_avail, 
         med_ppi, 
         by = c('aa_seq', 'deg1', 'deg2', 'pos1', 'pos2', 'max1', 'max2', 'max_pep', 'family'), 
         suffix = c('.avail', '.ppi'))

comp_signif$name <- 
gsub(comp_signif$family, pattern = 'PRM_0', replacement = '')



comp <- vector(mode='list', length = nrow(comp_signif))
for(i in 1:nrow(comp_signif)){
  
    x <- unlist(comp_signif[i, c('availability', 'binding')])
    comp[[i]] <- str_c( x, names(x), sep = ' ')
  
}

comp_signif$comp <- comp

#FigS7 <- 
ggplot(comp_signif[!is.na(comp_signif$binding) & !is.na(comp_signif$availability), ],
       aes(x = comp))+
  geom_bar(aes(fill = binding), color ='black', linewidth = 0.8, alpha = 0.85)+
  scale_y_continuous(transform = 'log2')+
  scale_x_upset(order_by = 'freq')+
  xlab('combinations')+
  scale_fill_manual(values=c( 'white','#92D74D','#3D3576'))+
  stat_count(
    geom = "text", colour = "black", size = 3,
    aes(x = comp, label = ..count..),
    position=position_stack(vjust=0.5))+
  t+
  theme(legend.position = 'inside', 
        legend.position.inside = c(0.75,0.8))+
  theme_combmatrix(combmatrix.label.text = element_text(size=8, color = 'black'),
                   combmatrix.label.extra_spacing = 5, 
                   combmatrix.panel.line.size = 1,
                   combmatrix.label.make_space = TRUE, 
                   combmatrix.panel.point.color.empty = 'grey', 
                   combmatrix.panel.point.alpha.empty =0.85, 
                   combmatrix.panel.point.size = 3)+
  guides(fill = guide_legend(title = 'Binding'))

saveRDS(FigS7, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/FigS7.rds')

Fig3F <- 
  ggplot(comp_signif[!is.na(comp_signif$binding) & !is.na(comp_signif$availability) & comp_signif$binding == 'stronger', ])+
  geom_bar(aes(x = availability, fill = availability), width = 0.6, color ='black', linewidth = 0.8)+
  scale_fill_manual(values=c('grey','#FE7F2D','#05004E'))+
  scale_x_discrete(limits = c('stronger', 'no significant difference', 'weaker'), 
                   labels = c('stronger', 'no significant \ndifference', 'weaker'))+
  t+theme(legend.position = 'none')

saveRDS(Fig3F, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/FigAvailF.rds')

# Fisher's exact test
# keep only sequences with the data for the 2 assays and remove stop codons

comp <- 
comp_signif[!(comp_signif$deg1 == '*' | comp_signif$deg2 == '*'), ]

comp <- 
  comp[!(is.na(comp$binding) | is.na(comp$availability)), ]

# Initialize variables
m <- sum(table(comp$binding)[c(2)])      # stronger PPI
n <- sum(table(comp$binding)[c(1,3)])        # weaker/similar PPI
k <-  sum(table(comp$availability)[c(2)])       # stronger avail
x <-  1:sum(table(comp[, c('binding', 'availability')])[c(2), c(2)]) # stronger PPI and stronger avail

# Use the dhyper built-in function for hypergeometric density
probabilities <- dhyper(x, m, n, k, log = FALSE)
probabilities

fisher.test(rbind(c(m,n),c(k,48)))

data <- data.frame( x = x, y = probabilities )
ggplot(data)+
  geom_density(stat="identity", aes(x, y))+
  geom_vline(xintercept = 49)



# Association looks different from the null hypothesis (which is that both variable are independent)

# save sequences with stronger or no difference in binding
comp%>%
  filter(binding != 'weaker')%>%
  select(aa_seq) -> specificity_pool

ref_sequence <- 
read_csv('~/PL_projects/sequencing/Screen_PPI_Avail_Feb2025/norm_sel_coef_ppi.csv')

specificity_pool <- 
ref_sequence[ref_sequence$aa_seq %in% specificity_pool$aa_seq, 'sequence']

dna_pool <- 
paste0('GGAGGTGGAGCTAGC', unlist(specificity_pool), 'aagcttattagttatgt')

pool <- 
tibble(pool_name = 'specificity_PBD', 
        sequence = dna_pool)

library(xlsx)
write.xlsx(pool, 'specificity_pool.xlsx')
