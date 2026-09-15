# AUTHOR : Pascale Lemieux
# analysis of sequence frequency changes 

#import packages and functions

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsci)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')

# set treshold of minimal reads at T0
min_reads = 20

#import datasets
reads_info <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/reads_info_both_PCA2.csv')
od_values <- read_csv('od_values.csv')
descript <- read_csv('sample_description.csv')
lib_descript <- read_csv('library_description.csv')

# compute number of generation per sample between timepoint with a starting OD od 0.15
od_values$gen_T1 <- 
unlist(lapply(od_values$OD_end_T1, get_generation, start_od = 0.15))

od_values$gen_T1_T2 <- 
  unlist(lapply(od_values$OD_end_T2, get_generation, start_od =0.15))

# sum the number of total generation over the assay
od_values$gen_T2 <- od_values$gen_T1 + od_values$gen_T1_T2

gen_values <- 
pivot_longer(od_values[,c(1,2,3,4,11, 13)], 
             names_to = 'timepoint', 
             cols = c(gen_T1, gen_T2))

gen_values$timepoint <- 
  gsub(pattern = 'gen_', replacement = '', gen_values$timepoint)

colnames(gen_values)[6] <- 'n_gen'

# replace na values for read count by 0
reads_info$number_both <- 
  replace_na(reads_info$number_both, 0)

# add 1 to the read count of each sequence
reads_info%<>%
  mutate(cor_count = number_both+1)

# long to wide format
seq_wide <- 
pivot_wider(reads_info, id_cols = c('sequence', 'aa_seq','family'), 
            names_from = c('timepoint', 'Pool', 'Replicate', 'Assay', 'condition'), 
            values_from = 'cor_count', values_fn = unique)


# select for sequences with > 20 reads at T0 in 2 out of 3 replicates
# keep only the info in the replicates with 20 reads or above at T0
for (pool in str_c('P', 1:6)) {
  
  sub_data <- seq_wide %>%
    select('sequence', 'aa_seq', 'family', grep(pool, colnames(seq_wide), value = TRUE))
  
  # Find the T0 columns (adjust pattern as needed)
  t0_cols <- grep('T0', colnames(sub_data), value = TRUE)
  
  T0 <- sub_data %>%
    filter(rowSums(across(all_of(t0_cols), ~ . > min_reads)) >= 2)
  
  # Keep only rows from original sub_data where sequence is in filtered T0
  sub_data <- sub_data[sub_data$sequence %in% T0$sequence, ]
  
  assign(pool, sub_data)
}

# compute fold change in each replicate considering only the sequence with reads>=20 at T0
pools <- list(P1, P2,P3, P4, P5, P6)

fold_change <- 
  lapply(pools, get_foldchange_per_condition, reads = min_reads)

names(fold_change) <- paste0('P', 1:6)

# pivot back to long with replicates
fold_change <- 
  lapply(fold_change, long_foldchange)

# associate generation number with log2 fold change
sample_gen <- 
merge(gen_values, 
  descript[, c(1,6)], 
  by = 'Sample')

sample_gen%<>%unique()

# compute selection coefficient
sel_coef <- 
  lapply(fold_change, get_selcoef, sample_gen)

sel_coef_df<-bind_rows(sel_coef)
# vizualise selection coefficient each sample computed using T1 or T2 as an end timepoint

ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/selection_coefficient_timepoint.png', width = 10, height = 5)
saveRDS(s_coef, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/selection_coefficient_timepoint.rds')

# write the selection coefficient for each pool in a dataframe format
write_csv(sel_coef_df, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/sel_coefficient_both_20.csv')



